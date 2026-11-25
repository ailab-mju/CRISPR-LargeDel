# CRISPR-LargeDel — 이식성/UB 수정 기록

작성일: 2026-08-17
목적: **g++-7 의존성을 제거하고, 컴파일러·최적화 레벨과 무관하게 동일한 결과가 나오도록** 만드는 것.
Rust/WASM 포팅의 정확성을 검증하려면 먼저 신뢰할 수 있는 기준 출력(golden fixture)이 필요합니다.

---

## 요약

| 항목 | 수정 전 | 수정 후 |
|---|---|---|
| g++ 13 빌드 | 실패 (3종 에러) | 성공 |
| clang 18 빌드 | 실패 | 성공 |
| 테스트 데이터 실행 (g++13 `-O3`) | **무한 루프** (10분+ 미완료, 8 GB 쓰레기 출력) | 0.12초 완료 |
| 테스트 데이터 실행 (g++13 `-O0`) | **SIGILL** | 정상 |
| ASan + UBSan | UB 4종 검출 | **경고 0** |
| `-O0/-O1/-O2/-O3` × g++/clang 출력 | 비교 불가 | **6개 빌드 전부 SAM 바이트 동일** |
| 원본 g++-7 결과 대비 `mut_only.sam` | — | **바이트 동일** |

핵심: **논문의 대형 결실 판정 결과(`*_mut_only.sam`)는 변하지 않았습니다.**

---

## 검증 결과

### 빌드 매트릭스

```
g++ 13.3    -O0 -O1 -O2 -O3   전부 빌드 성공
clang 18.1  -O0 -O2           전부 빌드 성공
```

6개 빌드 산출물의 SAM 파일 10종을 `cmp`로 비교 → **전부 바이트 단위 동일**.
보고 수치도 동일:

```
Control  sample Large Deletion Rate (without filtering) : 1.51%
Cas9 Mut sample Large Deletion Rate (removed false-pos) : 29.41%
Cas9 Mut sample small indel rate                        : 11.76%
```

### 새니타이저

```
g++ -O1 -fsanitize=address,undefined -fno-sanitize-recover=all
→ EXIT 0, stderr 0 bytes
```

### 원본(g++-7, Windows) 결과와의 대조

저장소에 커밋되어 있던 `result/HPRT1_TEST_11/` 를 기준으로:

| 파일 | 결과 |
|---|---|
| `HPRT1_TEST_mut_only.sam` | **바이트 동일** ← 대형 결실 판정 결과 |
| `HPRT1_TEST_site_ctr.sam` | **바이트 동일** |
| `HPRT1_TEST_site_mut_sm.sam` | **바이트 동일** |
| `HPRT1_TEST_ctr.sam` / `_mut.sam` | 정렬 좌표·CIGAR 기준 2,858건 일치 / 원본만 214건 / 신규만 330건 |
| `HPRT1_TEST_ctr_only.sam` | 원본 0바이트 → 수정 후 2,538바이트 |

**차이가 남는 부분에 대한 설명**

- `ctr.sam`/`mut.sam` 의 차이는 대부분 **`MISS:` 태그**입니다. `align_query()`의 `pre_chunk` 고속 경로는 `r.tags[1]`을 갱신하지 않아서, **직전 정렬 시도의 태그 값이 그대로 남습니다**(기존부터 있던 버그). 제어 흐름이 바뀌면서 남는 값이 달라진 것이고 정렬 자체는 동일합니다.
- `ctr_only.sam` 이 0 → 2,538 바이트가 된 것은, 원본에서는 쓰레기 반환값 때문에 control 쪽 split read가 통째로 누락됐기 때문입니다. 대조군 baseline이 새로 잡히는 것이므로 **§"영향 검토" 를 확인하세요.**

---

## 수정 내역 (7개 파일)

### 1. `align.h`, `post.h` — 빌드 실패 (이식성)

```cpp
+#include <cstdint>
+#include <array>
+#include <tuple>
```

최신 libstdc++는 `<cstdint>`를 더 이상 전이 include 하지 않아 `uint32_t` 미정의 에러가 납니다.

### 2. `Makefile` — 빌드 실패 (이식성)

```make
-CPP=/usr/bin/g++-7
+CXX ?= g++
+STD ?= -std=gnu++14
```

- 컴파일러 경로 하드코딩 제거
- **`-std=gnu++14` 고정이 필수**입니다. C++17부터 `std::data()`가 추가됐는데 모든 TU가 `using namespace std;`를 하고 있어서 사용자 정의 `struct data`와 충돌합니다 (`invalid types 'int[int]' for array subscript`). 근본 해결은 `struct data` 개명이지만, 기존 출력 보존을 위해 표준 버전 고정을 택했습니다.
- **주의**: `CPP`는 make의 내장 변수(C 전처리기)라 `CPP ?= g++`로 쓰면 `cc`로 남아 빌드가 깨집니다. 반드시 `CXX`를 쓰세요.
- `make clean` 타깃 추가

### 3. `func.cpp` — 반환값 누락 (**무한 루프의 주범**)

```cpp
 int long_inner_check(read &r1, read &r2, vector<read> *split, result &res)
 {
     ...
     // 마지막 3개 if 블록을 모두 통과하면 return 없이 함수 종료
+    return 0;
 }
```

호출부가 이 값으로 분기합니다:

```cpp
int res = long_inner_check(...);
if(res==-1){ return 0; }                    // 페어 폐기
if(res>=400){ align_query(...,1); ... }     // 강제 재정렬
```

즉 **레지스터 잔여값이 제어 흐름을 결정**했습니다. g++-7 `-O3`에서는 우연히 맞게 동작했고, g++13에서는 무한 루프, `-O0`에서는 SIGILL이 됐습니다.
이 경로는 "정상 페어, 대형 결실 없음"이므로 폐기도 재정렬도 트리거하면 안 됩니다 → `0`이 유일한 정답.

### 4. `align_wrap.cpp` — `read_fastq()` 버퍼 언더플로

`fgets()` 반환값을 무시하고 항상 `buf[strlen(buf)-1]`에 썼습니다. `strlen`이 0이면 **배열 앞 1바이트를 덮어씁니다**.

`strlen`이 0이 되는 경로가 실제로 존재합니다 — `char_to_num()`이 `'A'`를 바이트 `0x00`으로 인코딩하기 때문에, 이미 파싱된 리드 버퍼는 첫 염기가 A면 `strlen == 0`입니다. EOF에서 `fgets`가 버퍼를 그대로 두면 그 상태로 `buf[-1]`에 씁니다.

수정: 매 `fgets` 결과를 검사하고, 개행이 실제로 있을 때만 제거. 정상 FASTQ에 대해서는 동작이 동일합니다.

### 5. `align.cpp` — 경계 검사 4건

| 위치 | 문제 | 수정 |
|---|---|---|
| `side_align()` ×2 | `last`가 −1인 채로 정제 루프 진입 → `r.query[-1]`. 원래 `if(last==-1) return;` 가드가 **루프 뒤에** 있었음 | 가드를 루프 **앞**으로 |
| 좌측 확장 루프 | `st==0`에서 `g_ref_seq[-1]` | `&& lis_res[i].st>0` |
| 우측 확장 루프 | `ed+1`이 레퍼런스 끝을 넘어감 | `&& lis_res[i].ed+1<arg.ref_len` |
| microhomology `while(1)` | 하한 없이 4개 커서를 감소 → homopolymer/STR에서 음수 인덱스 | 경계에서 `break` |

### 6. `run_kmer_align.sh`

```sh
-mkdir $resdir
+mkdir -p $resdir
```

`result/` 가 없으면 전체 파이프라인이 실패했습니다.

---

## 영향 검토 — 논문 수치에 영향이 있나?

**대형 결실 판정 결과(`mut_only.sam`)는 바이트 단위로 동일합니다.** 이 파일이 보고된 Cas9 대형 결실률을 결정하므로, 논문의 핵심 수치는 그대로입니다.

다만 두 가지를 확인하시는 게 좋습니다:

1. **대조군 대형 결실(`ctr_only.sam`)이 0 → 8건**으로 바뀌었습니다. 원본에서는 UB 때문에 control split read가 통째로 누락된 것으로 보입니다. 대조군 baseline을 보고에 쓰셨다면 재확인이 필요합니다.
2. **전체 정렬률**: control 63.76%, mut 56.44%. 원본 대비 SAM 레코드가 3~5% 늘었습니다.

이 테스트 데이터는 2,500 read pair짜리 샘플이므로, **실제 논문 데이터셋으로 재실행해서 수치 변화 폭을 확인**하시길 권합니다.

---

## 아직 남겨둔 것 (의도적 — 알고리즘 동작을 바꾸므로)

수정하지 않았습니다. 판단이 필요합니다.

### A. `split[]` 이중 push 가능성

`align_pair()`에서:

```cpp
if(r1.map_f==1 && r2.map_f==1){
    int res = long_inner_check(...);   // 내부에서 map_f=2 설정 + split에 push
    if(res>=400){ align_query(...,1); align_query(...,1); }
}
if(r1.map_f==2 || r2.map_f==2){       // 위에서 2로 바뀌었으면 여기도 진입
    align_query(...,1); align_query(...,1);
    split[0].push_back(r1);            // 두 번째 push
    split[1].push_back(r2);
}
```

같은 페어가 `split[]`에 두 번 들어가 대형 결실이 중복 계수될 수 있습니다.

### B. `MISS:` 태그 staleness

`pre_chunk` 고속 경로가 `r.tags[1]`을 갱신하지 않아 이전 정렬 시도의 값이 남습니다. SAM 태그만의 문제이고 판정에는 영향 없습니다.

### C. `MIN_skip = 50` vs 논문의 ">100 bp"

코드는 50 bp 이상 gap을 split으로 인정합니다. 논문 본문 정의(>100 bp)와 다릅니다.

### D. `char t[200]` 초기화 누락

`make_ref()` / `char_to_num()`의 룩업 테이블이 A/C/G/T/N 인덱스만 설정되고 나머지는 초기화되지 않습니다. 레퍼런스나 리드에 소문자·IUPAC 모호 염기(R, Y, S, W …)가 있으면 쓰레기 값을 읽습니다. 현재 테스트 데이터에는 없지만 **실데이터에서는 터질 수 있습니다.**

### E. 고정 버퍼

`g_ref_seq[100000]`, `read.query[300]`, `mutual1[20000]`, `count1[20000]` — 전부 검사 없음.
레퍼런스 100 kb 초과, 리드 250 bp 초과, 클러스터 20,000개 초과 시 오버플로.

---

## 정정

이전 검토 문서에서 "k-means가 랜덤 초기화라 실행마다 결과가 달라질 수 있다"고 적었는데 **틀렸습니다.**
`dkm.hpp:79`에서 `rand_device()` 호출이 주석 처리되고 시드가 `22388974`로 하드코딩되어 있습니다.

```cpp
std::random_device rand_device;
std::linear_congruential_engine<uint64_t, ...> rand_engine( 22388974);
//		rand_device());
```

동일 바이너리 10회 반복 실행 결과 출력 해시가 전부 일치했습니다. **클러스터링은 결정론적입니다.**

---

## 재현 방법

```bash
make                       # g++ 기본
make CXX=clang++           # clang
make CXX=g++-7             # 원래 컴파일러도 그대로 동작

bash run_kmer_align.sh HPRT1 TEST 4652 4671 3 test_con test_mut 11

# 새니타이저 검증
g++ -O1 -g -std=gnu++14 -fsanitize=address,undefined -fno-sanitize-recover=all \
    -pthread -o kmer_align_san *.cpp
```

---

## 추가 수정 — 웹(WASM) 이식 과정에서 드러난 것

브라우저 빌드를 만들면서 네이티브에서는 우연히 가려져 있던 버그가 세 건 더 나왔습니다. 전부 본 소스에 반영했습니다.

### F. `printf` 인자 개수 불일치 (`align_wrap.cpp`)

```c
printf("Control sample Large Deletion Rate (wihtout filtering) : %.2lf%%\n",
       type.c_str(),                    // ← 포맷에는 %.2lf 하나뿐인데 인자가 둘
       (double)(...));
```

x86-64 SysV에서는 double이 `xmm0`, 포인터가 정수 레지스터로 따로 전달되어 **우연히 올바르게 출력**됐습니다. wasm은 가변 인자를 순서대로 스택에 쌓기 때문에 `printf`가 `char*`를 double로 읽어 `0.00`을 출력했습니다.

### G. `\%` 는 이스케이프가 아님 (`align_wrap.cpp`, 4곳)

```c
printf("... %.2lf\%\n", x);
```

`\%`는 C 이스케이프가 아니라서 `printf`는 `%\n`을 **변환 지정자로 해석**합니다(컴파일러도 `unknown conversion type character '\x0a'`로 경고). glibc는 그래도 `%`와 개행을 뱉었지만 musl(wasm)은 아무것도 출력하지 않아 로그 줄이 전부 붙어버렸습니다. `%%\n`으로 수정.

### H. k-means RNG가 표준 라이브러리에 의존 (`dkm.hpp`) — 가장 중요

두 가지가 겹쳐 있었습니다.

1. **엔진**: `std::linear_congruential_engine<uint64_t, 6364136223846793005, 1442695040888963407, UINT64_MAX>` 의 modulus가 2⁶⁴가 아니라 **2⁶⁴−1**입니다. `a*x`가 64비트를 넘치고, 이 나머지 연산을 libstdc++와 libc++가 다르게 처리합니다.
2. **분포**: `std::uniform_int_distribution` / `std::discrete_distribution`은 같은 엔진이라도 구현마다 다른 값을 냅니다(표준이 보장하지 않음).

결과적으로 **같은 시드인데 네이티브와 wasm이 서로 다른 초기 중심점**을 골랐고, 선택되는 클러스터 개수까지 달라졌습니다(같은 입력에 네이티브 9개 / wasm 6개).

수정: 재귀식을 2⁶⁴ 자연 wraparound로 직접 작성하고, 샘플링을 역-CDF 방식의 평범한 산술로 대체했습니다. 시드는 그대로입니다.

**영향**: 저장소 `test_data` 결과는 **변하지 않습니다**(`mut_only.sam` 바이트 동일, 29.41%). 시뮬레이션 데이터에서는 보고값이 최대 0.06 pp 움직였습니다(d20 19.92→19.98, d30 31.35→31.30). 정확도 오차 ±1.3 pp에 비하면 무시할 수준이고, 그 대가로 **어느 플랫폼에서든 같은 답**이 나옵니다.

### I. 편의 기능

- `argv[14]` 선택 인자 `full_sam`(기본 1). 0이면 전체 정렬 SAM 두 개를 생략합니다. 150k pair에서 파일당 약 110 MB라 브라우저에서는 실용적이지 않아 웹 빌드가 사용합니다. **기존 동작은 그대로**입니다.
- 위 옵션으로 출력 경로가 비었을 때를 위한 `fclose(NULL)` 가드.
- 단계별 `fflush(stdout)` — 브라우저에서 진행률 스트리밍용.

### 웹 빌드에서 확인된 것

네이티브 CLI와 브라우저 wasm이 **SAM 5개 파일 전부 바이트 동일**, 보고 수치도 전부 일치(29.41% / 11.76% / 56.44% / 63.76% / 1.51%). 속도는 wasm이 네이티브의 약 0.6배. 자세한 내용은 `web/README.md`.

> 주의: `exit()`이 대신 해주던 stdio flush가 wasm에는 없습니다(`EXIT_RUNTIME=0`). `align_wrap`이 대조군용으로 연 `f_out`/`f_out2`를 닫지 않고 처리군용으로 다시 열기 때문에, `fflush(NULL)` 없이는 그 파일들의 꼬리가 조용히 잘립니다. `web/wasm_entry.cpp`에서 처리합니다.

---

## 추가 (2026-09-20): usage 출력과 mutation pattern 요약

기존 동작·수치는 그대로이고, 출력만 추가됩니다. `test_data` 기준 **SAM 10종 전부 바이트 동일**, 기존 stdout 라인(`total_align`, `alignment_rate`, `Large Deletion Rate`, `small_indel rate`)도 문자열·순서 그대로라 `bench/evaluate.py` 와 `web/public/worker.js` 의 파싱은 영향받지 않습니다.

### 1. 인자 처리

이전에는 `Arg_Parse()` 가 `argv[1..12]` 를 검사 없이 읽어서, 인자 없이 실행하면 NULL 역참조로 바로 segfault 났습니다.

- `print_usage()` 추가. 인자 없음 / `-h` / `--help` → 사용법 출력 후 **exit 0**.
- 인자 12개 미만, 파일 열기 실패, 정수 파싱 실패, `target_end <= target_start`, `cleavage_offset` 범위 초과, `K` 가 5..31 밖, `out_dir` 부재·쓰기 불가 → 어느 번호 인자가 왜 틀렸는지 출력하고 사용법을 붙인 뒤 **exit 1**.
- `argc` 는 argv 의 NULL 종단을 세서 구하므로 `align_wrap(char*[])` 시그니처와 wasm 진입점은 그대로입니다. 결과는 `align_wrap::status` 로 전달되고 `main` / `run_analysis` 가 그대로 반환합니다.
- `run_kmer_align.sh` 도 인자 없이 실행하면 사용법을 출력합니다. SAM 정렬 단계는 빈 파일을 건너뛰도록 가드를 넣었습니다.

### 2. mutation pattern 별 count / %

카운터는 `struct result` 에 이미 쌓이고 있었지만 출력하는 `printf` 가 전부 주석 처리되어 있어서, 비율만 나오고 read count 는 SAM 줄 수를 직접 세야 했습니다.

- 실행 끝에 대조군·처리군 각각 `WT / Small insertion / Small deletion / Substitution / Large deletion` 의 read pair count 와 % 를 표로 출력합니다.
- 같은 내용을 `{gene}_{sample}_mutation_summary.tsv` 로도 씁니다.
- 분모는 기존 보고 수치가 쓰던 `tg_align_0 + 샘플 고유 대형 결실` 이라 5개 행이 100% 로 합해지고, 대조군과 공통인 클러스터(PCR 아티팩트)는 별도 행으로 뺍니다.

### 3. 부수 정리

- `post.cpp` 가 `extern int mut_only, ctr_only` 로, `align_wrap.cpp` 가 `long long int` 로 선언하던 타입 불일치(UB)를 `int` 로 통일했습니다.
- `write_read()` 가 only/common 카운트를 인자로 돌려주도록 바꿔, `ctr_only = mut_only` 로 값을 옮기던 호출 순서 의존을 없앴습니다. 값 자체는 동일합니다.
- `Makefile`: `-std=gnu++14` 고정(C++17 기본값에서 `struct data` 와 `std::data` 충돌), `g++-7` 없으면 `g++` 로 폴백, `make clean` 추가.

### 4. 절단부 끝단 soft clip 을 insertion 으로 계수

`make_cigar()` 는 정렬된 블록 두 개가 레퍼런스 좌표에서 겹칠 때만 `I` 를 만들기 때문에, 삽입이 read 끝에 걸리면 양옆 중 한쪽밖에 없어서 `I` 가 아니라 soft clip 이 됩니다. 그 결과 read 끝 삽입은 전부 WT 로 집계되고 있었습니다.

- `read_classify()` 에 `classify_clip()` 추가: **5 bp 이상**의 끝단 soft clip 경계가 절단부 **±13 bp** 안이면 `r.sclip` 을 세웁니다 (임계값은 `align.cpp` 의 `MIN_CLIP` / `CLIP_WIN`).
- `long_inner_check()` 가 `tg_clip` 으로 세고 해당 read 를 `_site_{ctr,mut}_clip.sam` 에 씁니다. 판정 우선순위는 `I` > `D` > clip > 치환입니다.
- SAM 의 첫 optional field(원래 항상 비어 있던 자리)에 `CLIP:<길이>:<L|R>` 태그를 붙입니다.
- 요약표·TSV 에 `Terminal insertion (clipped)` 행 추가. **기존 `small_indel rate` 계산식에는 넣지 않았습니다** — 이전에 보고한 수치가 그대로 유지됩니다.

시뮬레이션 검증 (HPRT1, 절단부 4655, 100 pair): 절단부에서 끝나는 8/10/20/30 bp 끝단 삽입 → 전부 `Terminal insertion (clipped)` 100%. 절단부에서 60 bp 떨어진 동일한 clip → WT 유지. 중간 삽입 10 bp → 기존대로 `Small insertion`.

한계 두 가지:
- 정렬 확장이 mismatch 1개를 허용하며 밀고 들어가므로 N bp 끝단 삽입은 대략 `(N-1)S`~`(N-2)S` 로 나타납니다. 6 bp 미만 삽입은 chunk 병합 과정에서 짧은 `D` 로 흡수되어 `Small deletion` 으로 가는 경우가 많습니다.
- read 길이의 20% 를 넘는 끝단 clip 은 `make_cigar()` 의 `r_match < 0.8*query_len` 에 걸려 pair 가 통째로 버려집니다 (150 bp read 기준 32 bp 이상).

`test_data` 회귀: SAM 10종 중 9종 바이트 동일, `_mut.sam` 만 1 pair 에 `CLIP:8:L` 태그가 추가되어 다릅니다(그 pair 는 절단부를 가로지르지 않아 어떤 카테고리에도 들어가지 않습니다). 보고 수치 4종(29.41% / 11.76% / 56.44% / 63.76%)은 전부 동일합니다.

### 5. 단일 분모 · D 앵커 최소 길이 · read 이름에 판정 기록

실데이터(`Issue/PRR5/res_v2`)에서 "100% 편집된 라인이 WT 47% / 결실 51% 로 나온다"는 문제를 추적하다 나온 세 가지 수정입니다.

**(a) 판정 분모를 pair 바깥 좌표 → 실제 커버리지로**

`long_inner_check()` 의 절단부 통과 판정이 pair 의 **바깥 좌표**(`l_p`..`r_p`)를 쓰고 있었습니다. mate 두 개가 절단부를 사이에 두고 양옆에 앉으면 **둘 다 절단부를 전혀 읽지 않았는데도** 조건을 통과하고, 양쪽 모두 완벽 정렬되니 WT 로 계수됐습니다. PRR5 처리군에서는 WT 로 분류된 10,383 pair 중 접합부를 실제로 덮는 게 **1개**였습니다.

이제 **최소 한쪽 mate 가 절단부 ±13 을 실제로 덮거나(`r.sclip` 포함)** 해야 계수합니다. 걸쳐만 있고 못 덮은 pair 는 `tg_nocov` 로 따로 세어 `Not covered (excluded)` 행에 표시하고 모든 비율에서 뺍니다. 대형 결실은 이 검사 이전 분기에서 빠져나가므로 영향 없습니다.

분모가 하나로 통일되면서 기존 `Large Deletion Rate` / `small_indel rate` 두 줄도 같은 분모를 쓰게 됩니다. **이전에 보고된 수치와 달라집니다.**

**(b) `D` 는 양옆에 `MIN_match` 이상 블록이 있을 때만**

`side_align()` 이 첫/마지막 블록을 최대 30 bp 옮겨 mismatch 를 줄이는데, 길이 검사 없이 결과를 push 합니다. 그래서 절단부 하류에 우연히 2~3 bp 만 맞아도 큰 `D` 가 만들어집니다. PRR5 에서는 편집 allele 의 새 접합부가 우연히 26 bp 하류 레퍼런스와 11 bp 일치해서, 삽입이 `139M26D3M9S` 같은 "26 bp 결실"로 둔갑했습니다 (26D read 9,463개 중 7,998개가 `3M` 앵커).

`prune_short_ends()` 로 양 끝의 `MIN_match`(13 bp) 미만 블록을 버립니다. `chunk_align()` 이 이미 시드 청크에 같은 하한을 걸고 있으므로 `side_align()` 산출물만 영향받습니다. 남은 부분은 soft clip 이 되고 `classify_clip()` 이 끝단 삽입으로 잡습니다.

**(c) read 이름에 판정 결과 기록**

`print_read()` 가 QNAME 뒤에 `|<패턴>` 을 붙입니다: `WT` `INS` `DEL` `SUB` `CLIP` `LDEL_ONLY` `LDEL_COMMON`, 그리고 전체 정렬 SAM 에만 나오는 `NOCOV`(걸쳤지만 못 덮음) `OUT`(절단부와 무관). 두 mate 는 항상 같은 접미사를 받으므로 pairing 에 영향 없습니다.

**검증 — `Simulation_data/subset/noPCR` (HPRT1, ART 시뮬, read 이름에 정답 allele)**

정답 클래스별 판정 분포 (판정된 pair 기준):

pair 단위 정답 클래스 x 판정 (판정된 pair 기준):

| 정답 | 수정 전 | 수정 후 |
| --- | --- | --- |
| DEL 1~10 bp | WT **8.5~20.4%** / DEL 75~88% | WT **0.0%** / DEL 80~88% / CLIP 8~15% |
| DEL 100 bp | WT 15.2% / LDEL_ONLY 79.0% | WT **0.0%** / LDEL_ONLY **93.2%** |
| DEL 500~3000 bp | LDEL_ONLY 96~99% | LDEL_ONLY 96~99% (변화 없음) |
| INS 1~5 bp | WT 12.2~17.2% / INS 71~81% | WT **0.0%** / INS 79~87% / CLIP 12~19% |
| INS 10 bp | WT 31.7% / INS 48.3% | WT **0.0%** / INS 48.2% + CLIP 48.2% |
| INS 100 bp | INS **0%** (DEL 27.3%, SUB 27.3%, LDEL_COMMON 27.3% 로 흩어짐) | **CLIP 87.1%** |
| WT | WT 90.8% / SUB 4.9% | WT 89.7% / SUB 5.6% |

즉 **모든 변이 클래스에서 WT 오분류가 0 이 되었고**, 기존에 아예 안 잡히던 100 bp 끝단 삽입이 87% 잡힙니다. 진짜 WT 의 정분류율은 그대로입니다.

`withPCR` 서브셋(키메라 포함)에서도 같은 경향이며, INS 100 bp 는 CLIP 96.2% 입니다.

남은 한계 — 시뮬레이션에서 **확인된** 것:
- `CLIP` 은 이름은 "Terminal insertion" 이지만 **끝단 결실도 같이 잡습니다**. `_site_mut_clip.sam` 안에 정답이 DEL 1~10 bp 인 read 가 126 개 들어 있고 CIGAR 는 전부 `8S142M` `141M9S` 같은 끝단 clip 입니다. read 끝에서는 삽입과 결실이 똑같이 soft clip 으로 보이므로 원리상 구분이 안 됩니다.
- 5000 bp 이상 대형 결실은 only/common 분할이 무너집니다 (DEL 5000: COMMON 40~52%, DEL 7000: COMMON 100%). 다만 판정된 pair 가 21~35 개로 적습니다.
- 진짜 WT 의 1.7~2.4% 가 `LDEL_ONLY` 로 잡힙니다 (대형 결실 위양성 바닥값).

시뮬레이션에서 **재현되지 않은** 것:
- PRR5 실데이터의 "처리군 대형 결실 대부분이 `common` 으로 제외됨" (16,668 / 18,462 = 90%) 은 시뮬에서 나타나지 않습니다. noPCR 9.4%, withPCR 10.5% 입니다. 시뮬의 키메라는 두 시료에 동일 좌표로 심겨 있어 k-means 가 의도대로 제거하는 반면, PRR5 의 거동은 별개 현상이며 아직 미해명입니다.

### 6. SAM 헤더 · 유효한 태그 · htslib 링크로 BAM 직접 출력

지금까지 출력 SAM 에는 `@HD`/`@SQ` 헤더가 없었고, QUAL 뒤에 빈 칼럼이 하나 붙었으며(`tags[0]` 이 항상 빈 문자열), 태그가 `MISS:2:0` `Dist:50000` `D_size:151` 처럼 SAM 규격(2글자 태그 + 타입 문자)을 어겼습니다. 그래서 samtools 로 바로 읽을 수 없어 `Issue/fix_sam.sh` 로 후처리해야 했습니다.

**출력 계층을 새로 만들었습니다 (`samout.h` / `samout.cpp`).**

- htslib 을 링크해 `kmer_align` 이 **좌표 정렬 BAM + `.bai` 를 직접** 씁니다. 파이프라인 어디에도 samtools 호출이 없고, 수 GB 짜리 중간 SAM 도 만들지 않습니다.
- 정렬은 **위치 버킷 방식**입니다. 레코드를 256개 버킷(각 64 KB 버퍼)으로 나눠 임시 파일에 흘려보내고, close 시 버킷 순서대로 읽어 정렬해 내보냅니다. 깊이와 무관하게 메모리가 한정됩니다.
- `sam_parse1()` 로 텍스트 레코드를 `bam1_t` 로 변환하므로 CIGAR·SEQ·태그 인코딩을 손으로 짤 필요가 없고, `sam_idx_init()`/`sam_idx_save()` 가 쓰는 동안 인덱스를 만듭니다.
- htslib 이 없으면 같은 코드가 정렬만 하고 SAM 을 씁니다 — 웹(wasm) 빌드가 이 경로를 씁니다(`viewer.js` 가 SAM 텍스트를 읽으므로).
- 15번째 인자 `out_fmt` = `bam`(기본) / `sam` / `both`. 래퍼는 `OUT_FMT=` 환경변수로 넘깁니다.
- Makefile 이 htslib 을 자동 탐지하고 어느 쪽인지 출력합니다. `HTS=0` 으로 끄고, `HTS_CFLAGS`/`HTS_LDFLAGS`/`HTS_LIBS` 로 경로를 지정합니다.

**헤더와 태그**

- 모든 출력에 `@HD VN:1.6 SO:coordinate`, `@SQ SN:<gene> LN:<ref_len>`, `@PG ID:kmer_align CL:<명령줄>`.
- `print_read()` 를 `sam_record()` 로 바꾸고 빈 태그 슬롯을 건너뜁니다 — 규격을 깨던 `\t\t` 가 사라집니다.
- 태그: `XM:i:`(총 mismatch) `XT:i:`(절단부 ±5 내) `XD:i:`(절단부까지 거리) `XL:i:`(정렬 reference 스팬 = split read 의 결실 크기) `XC:i:`(끝단 clip 길이) `XS:A:`(clip 방향 L/R) `XV:Z:`(판정 패턴) `YC:Z:`(클러스터 색, 기존 유지).

**함정 하나** — `sam_idx_init(fp,hdr,0,(base+".bam.bai").c_str())` 처럼 임시 문자열의 `c_str()` 을 넘기면 htslib 이 그 포인터를 보관하는데 임시가 즉시 소멸해 **인덱스가 조용히 안 만들어집니다**(반환값은 0). 경로 문자열을 지역 변수로 잡아둬야 합니다.

**부수 수정**

- `web/public/viewer.js` 의 `parseTags()` 가 `TAG:TYPE:VALUE` 의 타입 문자를 안 버리고 있었습니다(`XM` 이 `"i:2"` 로). 고쳤고 `app.js` 툴팁에 판정·끝단 clip 을 추가했습니다.
- `web/build.sh` 에 `samout.cpp` 추가, `wasm_entry.cpp` 가 `samout_set_format(OUT_SAM)` 을 호출합니다.
- `run_kmer_align.sh` 에서 `sort -k4n` 12줄과 samtools 단계를 전부 걷어냈습니다. 바이너리가 정렬해서 냅니다.
- `bench/evaluate.py`, `Simulation_data/scripts/s7_score_calls.py` 가 `@` 헤더 줄을 건너뛰고, BAM 만 있으면 `samtools view -h` 로 읽습니다.
- `Issue/fix_sam.sh` 는 구버전 결과(`res_v2` 이전) 전용으로 남겼습니다.

**`-1*` CIGAR (BAM 전환하고 나서야 드러난 버그)**

`make_cigar()` 는 "이 read 는 정렬 안 됨"을 `{len:-1, type:'*'}` 한 칸으로 표시합니다. 이걸 그대로 찍으면 CIGAR 가 `-1*` 이 되는데, SAM 규격상 불법입니다. 텍스트 SAM 시절에는 아무도 파싱을 안 해서 조용히 넘어갔지만 htslib 은 레코드째로 거부합니다(`[E::parse_cigar] CIGAR length invalid`).

어떻게 새어 나왔냐면 — `align_pair()` 가 `map_f==2` 경로에서 재정렬을 돌리는데, 재정렬이 한쪽 mate 를 무효화해도 여전히 `1` 을 반환합니다. 그래서 `process()` 가 그 pair 를 전체 정렬 출력에 씁니다. POS/PNEXT/TLEN 도 직전 정렬 시도에서 남은 값입니다.

- `read_aligned()` 추가 — `map_f==-1` 이거나 CIGAR 가 sentinel 이면 false.
- `process()` 가 두 mate 중 하나라도 정렬 실패면 **pair 전체를 전체 정렬 출력에서 제외**하고, 그 수를 `incomplete pairs not written to the whole-alignment output:` 로 보고합니다. 이전에는 htslib 이 나쁜 쪽만 버려서 **짝 없는 mate 가 BAM 에 남았습니다**.
- 방어적으로 `sam_record()` 도 sentinel 을 `*` 로 씁니다.
- `total_align` / `alignment_rate` 는 건드리지 않았습니다 — 보고 수치 불변.

`stress_large` 로 검증: 경고 8건 → 0건, 제외된 pair 6(대조군)/8(처리군), 보고 수치 전부 동일, 전체 정렬 BAM 레코드만 20개 감소(= 14 pair × 2 − htslib 이 이미 버리던 8개).

**`make_ref()` 의 `fscanf` 반환값** — 경고로만 보이던 것이 실제 버그였습니다. FASTA 가 비었거나 서열 줄이 없으면 `g_ref_seq` 에 이전 내용이 남은 채로 실행이 계속됐습니다. 이제 에러를 내고 종료합니다.

**검증** — `test_data` 12종 전부 `samtools quickcheck` 통과, `.bai` 생성, 레코드 수가 이전 SAM 빌드와 일치(3188/2822/48/20/…), 좌표순 정렬 확인, `both` 모드에서 BAM 과 SAM 내용 일치, htslib 없는 빌드는 정렬된 SAM 출력, 보고 수치 변동 없음(27.78 % / 11.11 %).

---

## 2026-09-21 — large deletion 필터: 논문 대조 + 시뮬레이션/PRR5 실측

### 먼저, 논문에 400 bp 는 없습니다

Hwang & Lee et al., *Nat. Biomed. Eng.* **9**:79–92 (2025) 기준:

- p.80 — "We defined a small insertion and a small deletion event as shorter than 100 bp and a **large-deletion event as larger than 100 bp**."
- p.89 — "Read pairs with a **deletion length of 100 bp or more** … are determined as large deletions."
- p.89 — read pair 3분류: "(1) all read pairs are skewed to the cleavage site, **(2) all read pairs are mapped without splitting and the cleavage site passes between reads**, and (3) one of the reads is split and passes through the cleavage site. … Read pairs in categories (2) and (3) are considered as potential CRISPR-derived variant candidates."
- p.89 — 분모 N_t = "the total number of reads passing the cut site (±13 base pairs), **including large deletions**" → 이미 반영된 단일 분모 방식과 일치합니다.

즉 **pair 간격 기반 호출(category 2)은 논문 알고리즘의 절반이라 뺄 수 없고**, `400` 은 논문에 근거가 없는 구현상의 숫자입니다.

### 그런데 400 을 100 으로 낮추면 안 됩니다 (실측)

논문의 100 bp 는 **deletion 길이** 정의입니다. split 되지 않은 pair 가 알려주는 건 **두 mate 사이 거리(insert size)** 뿐이고, ~300 bp 라이브러리에서는 정상 WT fragment 도 수백 bp 간격이 납니다. 그래서 `400` 은 실은 *deletion 길이 임계*가 아니라 **insert size 하한**입니다. 낮춰서 측정해보면:

| | small DEL 민감도 | large DEL 민감도 | WT→LDEL 오분류 | PRR5 LDEL% |
|---|---|---|---|---|
| min_gap 400 (기존) | 69.29 % | 51.27 % | 0.14 % | 8.79 % |
| min_gap 100 (논문 숫자) | 59.39 % | 42.95 % | **0.61 %** | **40.14 %** |

PRR5 는 large deletion 이 원래 없는 line 이므로 40.14 % 는 전부 오검출입니다. **min_gap 은 400 유지**하고, 의미를 주석·usage 에 명시한 뒤 인자로 노출했습니다.

같은 이유로 split read 쪽 `MIN_skip` 도 50 유지입니다. 100 으로 올리면 50–100 bp deletion 민감도가 69.29 % → 40.34 % 로 반토막 납니다.

### 채택한 필터 2개

**① cut_between (기본 켬)** — 논문 category (2) 문구 그대로, cleavage site 가 두 mate **사이**에 들어와야 합니다. 기존 코드는 gap 이 protospacer 와 겹치기만 하면 통과시켰습니다. `cleavage_offset` 이 크면(PRR5 는 16) cut 이 gap 밖에 있는데도 통과합니다.

**② anchor_min (기본 13)** — split read large deletion 의 anchor 품질. 정렬 블록마다 `길이 − 4 × mismatch` 를 계산해 **가장 나쁜 블록**의 점수를 씁니다. 긴 본체 블록은 수백 점이라 사실상 anchor 검사입니다. 13 이면 13 bp anchor 는 완전 일치해야 하고 17 bp anchor 는 mismatch 1개까지 허용됩니다.

`side_align()` 은 끝 블록을 30 bp 안에서 mismatch 가 줄어드는 아무 위치로나 옮기고, `prune_short_ends()` 는 13 bp 만 넘으면 살립니다. PRR5 의 유령 `26D` 가 정확히 이렇게 만들어졌습니다. 점수는 항상 `XA:i:` 태그로 나가므로 사후 필터링도 가능합니다.

### 결과

**시뮬레이션** (`stress_large`, HPRT1, 표본당 77k pair, read 이름의 ground truth 로 채점. 분모 = target window 에 닿은 pair)

| 설정 | small DEL 50–300 | large DEL 500+ | WT→LDEL | CHI→only | CHI→common |
|---|---|---|---|---|---|
| 기존 | 69.29 % | 51.27 % | 0.14 % | 43.32 % | 5.26 % |
| + cut_between | 69.29 % | 51.23 % | 0.14 % | 43.32 % | 5.26 % |
| + anchor 13 | 68.98 % | 52.26 % | 0.11 % | 43.32 % | 5.06 % |
| **둘 다 (새 기본값)** | **68.98 %** | **52.23 %** | **0.11 %** | 43.32 % | 5.06 % |
| + anchor 16 | 67.87 % | 51.50 % | 0.11 % | 42.71 % | 5.06 % |

민감도 손실 0.31 pp, large deletion 은 오히려 0.96 pp 상승했습니다.

**PRR5** (실데이터 절반: 처리군 1.27 M pair, 대조군 0.88 M pair. 이 line 은 +16 bp insertion 하나뿐이고 large deletion 이 없으므로 LDEL 은 전부 artefact)

| 설정 | analyzed | LDEL | LDEL % | 말단 insertion | 그 비율 | mut_only | mut_common | ctr LDEL % |
|---|---|---|---|---|---|---|---|---|
| 기존 | 5975 | 525 | 8.79 % | 5432 | 90.91 % | 525 | 8688 | 6.44 % |
| + cut_between | 5691 | 241 | 4.23 % | 5432 | 95.45 % | 241 | 8909 | 5.20 % |
| + anchor 13 | 5765 | 315 | 5.46 % | 5432 | 94.22 % | 315 | 742 | 3.56 % |
| **둘 다 (새 기본값)** | 5767 | **317** | **5.50 %** | **5432** | **94.19 %** | 317 | 677 | 2.50 % |
| + anchor 16 | 5618 | 168 | 2.99 % | 5432 | 96.69 % | 168 | 579 | 0.41 % |

**말단 insertion 이 모든 설정에서 5432 로 동일합니다** — 진짜 변이는 어느 필터도 건드리지 않습니다. 오검출만 525 → 317 (−40 %) 로 줄고, 그만큼 진짜 변이 비율이 90.91 % → 94.19 % 로 올라갑니다.

**PRR5 mut_only 병리 프로파일** (지적하신 4가지가 실제로 줄었는지)

| | 기존 | cut_between | anchor 13 | 둘 다 |
|---|---|---|---|---|
| pair 수 | 525 | 241 | 315 | 317 |
| 양쪽 다 D≥50bp 없음 (순수 간격 호출) | 389 (74.1 %) | 111 (46.1 %) | 241 (76.5 %) | 231 (72.9 %) |
| mate 가 삭제 구간 안에 있음 | 3 | 3 | 1 | 1 |
| RF 방향 | 84 (16.0 %) | 78 | 56 | 70 |
| anchor < 20 bp 인 split read | 41/136 (30.1 %) | 38/130 | **4/74 (5.4 %)** | 6/86 (7.0 %) |
| anchor 중앙값 | 28 bp | 28 bp | **44 bp** | 44 bp |

두 필터가 서로 다른 고장을 잡습니다. cut_between 은 간격 호출(389 → 111), anchor_min 은 가짜 split(41 → 4).

### 기각한 필터 2개 — 오히려 나빠집니다

**RF 방향 제거** 와 **모순 pair 제거**(mate 가 삭제 구간 안) 는 둘 다 artefact 를 **대조군에서도 똑같이** 지웁니다. 그런데 논문 3단계(k-means 로 양쪽 공통 cluster 제거)는 **대조군이 그 artefact 를 그대로 들고 있어야** 성립합니다. 지워버리면 원래 `common` 으로 가던 처리군 pair 가 `only` 로 넘어옵니다.

| | 시뮬 CHI→common | 시뮬 CHI→only | PRR5 mut_only |
|---|---|---|---|
| 기존 | 5.26 % | 43.32 % | 525 |
| 모순 pair 제거 | 2.83 % | 45.55 % | **684** |
| RF 제거 | **0.00 %** | 47.17 % | **3552** |

RF 방향 자체는 여전히 chimera 신호가 맞습니다(mut_only 16.0 % vs mut_common 60.4 %). 다만 **버리는 게 아니라 대조군에 남겨 상쇄시키는 게** 이 알고리즘의 작동 방식입니다.

### 코드 변경

| 파일 | 변경 |
|---|---|
| `arg_parse.h` | `Args` 에 `min_gap`, `anchor_min`, `cut_between` 추가 |
| `arg_parse.cpp` | 선택 인자 16/17/18 파싱 + 검증, usage 에 설명 추가 (400 이 insert size 하한이라는 경고 포함) |
| `func.cpp` | `long_inner_check()`: 하드코딩 `400` → `arg.min_gap`, cut_between 게이트 추가, 근거 주석 |
| `align.cpp` | `read_classify()`: anchor 점수 계산 + `XA:i:` 태그 + `arg.anchor_min` 필터 |
| `align_wrap.cpp` | tag 슬롯 4 → 5, `res>=400` → `res>=arg.min_gap` |
| `run_kmer_align.sh` | `MIN_GAP` / `ANCHOR_MIN` / `CUT_BETWEEN` 환경변수 전달 + 도움말 |

기본값은 `min_gap=400`, `anchor_min=13`, `cut_between=1`. 인자를 안 주면 위 표의 "둘 다" 열이 나옵니다. 예전 동작이 필요하면 `ANCHOR_MIN=0 CUT_BETWEEN=0`.

### 검증

- `test_data` 소량 세트: BAM + `.bai` 정상 생성, `both` 모드 정상, 크래시 없음.
- 시뮬레이션·PRR5 모두 **프로덕션 빌드 결과가 실험 빌드(CBA13)와 바이트 단위로 동일**.
- **주의**: 부분 재빌드 시 `-DUSE_HTSLIB` 유무가 섞인 `.o` 가 남으면 `munmap_chunk(): invalid pointer` 로 죽습니다. htslib 설정을 바꿀 때는 `make clean` 후 전체 재빌드하세요.
