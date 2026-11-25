/* WebAssembly entry point for the ExCas-Analyzer core (kmer_align).
 *
 * The native tool is driven entirely through argv, so this just rebuilds an
 * argv array from typed arguments and hands it to align_wrap - the analysis
 * code itself is untouched and behaves exactly as the native build.
 *
 * A fresh Module instance is created per run on the JS side, because the
 * core keeps its state in file-scope globals (g_ref_seq, arg, mut_only,
 * ctr_only, skip_read ...) and never resets them.
 *
 * Based on /home/seb/web/wasm_entry.cpp (2026-10-02), extended with the two
 * arguments docs/worker.js passes since v20261023: clip_rescue and
 * common_mode.  Options after #24 (the detection / speed switches of the
 * cli_detection and cli_detect_speed builds) are not passed, so they keep
 * their defaults, which are on.
 */
#include <emscripten.h>
#include <stdio.h>
#include <string>
#include <vector>
#include "align_wrap.h"
#include "../samout.h"

extern "C" {

/* Positional argv layout (see Arg_Parse / print_usage):
 *   1-5  ref ctr1 ctr2 mut1 mut2      6-8  target_start target_end cleavage
 *   9-10 gene sample                  11   K        12 outdir
 *   13   dist_cut                     14   full_sam
 *   15   out_fmt (always "sam" here - docs/viewer.js reads SAM text)
 *   16   min_gap   17 anchor_min   18 cut_between
 *   19   ins_fix   20 ins_max      21 report
 *   22   large_del_len   23 clip_rescue   24 common_mode (0 = local, 1 = kmeans)
 * Every argument is passed explicitly, so the page's defaults must mirror the
 * CLI's for the two to agree. */
EMSCRIPTEN_KEEPALIVE
int run_analysis(const char* ref,
                 const char* ctr1, const char* ctr2,
                 const char* mut1, const char* mut2,
                 int target_start, int target_end, int cleavage_offset,
                 const char* gene, const char* sample,
                 int k, const char* outdir,
                 int dist_cut, int full_sam,
                 int min_gap, int anchor_min, int cut_between,
                 int ins_fix, int ins_max, int report,
                 int large_del_len, int clip_rescue, int common_mode)
{
    std::vector<std::string> nums = {
        std::to_string(target_start), std::to_string(target_end),
        std::to_string(cleavage_offset), std::to_string(k),
        std::to_string(dist_cut), std::to_string(full_sam ? 1 : 0),
        std::to_string(min_gap), std::to_string(anchor_min),
        std::to_string(cut_between), std::to_string(ins_fix),
        std::to_string(ins_max), std::to_string(report),
        std::to_string(large_del_len), std::to_string(clip_rescue ? 1 : 0),
    };
    const std::string fmt = "sam";
    const std::string mode = common_mode ? "kmeans" : "local";

    std::vector<char*> argv;
    auto push = [&argv](const char* v){ argv.push_back(const_cast<char*>(v)); };
    push("kmer_align");
    push(ref); push(ctr1); push(ctr2); push(mut1); push(mut2);
    push(nums[0].c_str()); push(nums[1].c_str()); push(nums[2].c_str());
    push(gene); push(sample); push(nums[3].c_str()); push(outdir);
    push(nums[4].c_str()); push(nums[5].c_str());
    push(fmt.c_str());
    for (size_t i = 6; i < nums.size(); i++) push(nums[i].c_str());
    push(mode.c_str());
    argv.push_back(nullptr);

    samout_set_format(OUT_SAM);   /* the page reads SAM text out of the FS */
    align_wrap run(argv.data());
    /* Arg_Parse validates the arguments now and prints a usage block on bad
       input; run.status is non-zero in that case and nothing was analysed. */

    /* The native program relies on exit() to flush and close its output
       streams.  The wasm runtime is kept alive (EXIT_RUNTIME=0), so flush
       every open output stream here. */
    fflush(NULL);
    return run.status;
}

}  /* extern "C" */

/* emscripten still wants a main; the runtime is started with INVOKE_RUN=0 so
   this never executes. */
int main() { return 0; }
