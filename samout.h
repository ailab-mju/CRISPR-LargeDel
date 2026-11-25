#ifndef SAMOUT_H
#define SAMOUT_H
#include <string>

/* Coordinate-sorted, indexed alignment output, done inside kmer_align.
 *
 * Records arrive in FASTQ order but IGV needs them in reference order with a
 * .bai, so they are bucketed by position into small temp files while the run
 * proceeds and merged back in order at close().  Memory stays bounded (one
 * 64 KB buffer per bucket) no matter how large the run is, and the whole-
 * alignment output never has to exist as a multi-GB text file.
 *
 * Built WITHOUT htslib (the wasm build) this still sorts, but only emits SAM.
 */
struct read;

#define OUT_BAM  0
#define OUT_SAM  1
#define OUT_BOTH 2

struct SamOut;

void    samout_set_format(int fmt);      /* OUT_BAM by default */
int     samout_format(void);
bool    samout_have_bam(void);           /* false when built without htslib */

/* `base` is the output path WITHOUT extension. Returns NULL when base is
   empty, so callers keep using the same `if (out != NULL)` guard. */
SamOut *samout_open(const std::string &base);
void    samout_write(SamOut *o, read &r);
void    samout_close(SamOut *o);         /* sorts, writes .bam/.sam + .bai */

#endif
