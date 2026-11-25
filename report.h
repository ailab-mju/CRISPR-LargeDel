#ifndef KA_REPORT_H
#define KA_REPORT_H
#include <string>

/* HTML 리포트.

   실행 중에는 samout_write() 를 지나가는 레코드를 보고 인델 길이, 대형 결실
   길이, 8 bp 이상 삽입 서열을 집계합니다.  분류 코드는 건드리지 않습니다 -
   파일로 나가는 것과 정확히 같은 레코드를 세므로 리포트와 SAM 이 어긋날 수
   없습니다.

   실행이 끝나면 report_finish() 가
     <outdir>/<gene>_<sample>_report.html   이 실행 하나의 리포트
     <outdir>/<gene>_<sample>_report.tsv    합치기용 집계
   를 씁니다.  여러 실행을 합치거나 수정 전후를 비교하려면

     kmer_align report -o all.html [--compare DIR...] -- DIR...

   가 각 폴더의 _report.tsv 를 읽어 같은 템플릿으로 하나로 묶습니다. */

enum { RR_NONE=0, RR_MUT_SITE, RR_CTR_SITE, RR_MUT_LDEL, RR_CTR_LDEL };

int  report_role_for(const std::string &base);             /* samout_open 에서 */
void report_observe(int role, const std::string &sam_line); /* samout_write 에서 */

struct ka_pat { int wt,ins,del,sub,clip,ldel,lcommon,nocov,total,aligned; };
void report_finish(const ka_pat &ctr, const ka_pat &mut);

int  report_main(int argc, char **argv);   /* `kmer_align report ...` */
#endif
