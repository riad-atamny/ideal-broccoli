// File: stdinc.h
// David E. Taylor
// Applied Research Laboratory
// Department of Computer Science and Engineering
// Washington University in Saint Louis
// det3@arl.wustl.edu
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "uint128_t.h"

typedef char bit;
// const int false = 0;
// const int true = 1;

const int Null = 0;
const int BIGINT = 0x7fffffff;
const int EOS = '\0';

inline int max(int x, int y) { return x > y ? x : y; }
inline double max(double x, double y) { return x > y ? x : y; }
inline int min(int x, int y) { return x < y ? x : y; }
inline double min(double x, double y) { return x < y ? x : y; }
inline int abs(int x) { return x < 0 ? -x : x; }
// inline bit isdigit(int c) { return (c >= '0') && (c <= '9'); }

inline void warning(char* p) { fprintf(stderr,"Warning:%s \n",p); }
inline void fatal(char* string) {fprintf(stderr,"Fatal:%s\n",string); exit(1); }

double pow(double,double);
double log(double);

// long random(); double exp(double),log(double);

// Return a random number in [0,1] 
inline double randfrac() { return ((double) random())/BIGINT; }

// Return a random integer in the range [lo,hi].
// Not very good if range is larger than 10**7.
inline int randint(int lo, int hi) { return lo + (random() % (hi + 1 - lo)); }

// Return a random number from an exponential distribution with mean mu 
inline double randexp(double mu) { return -mu*log(randfrac()); }

// Return a random number from a geometric distribution with mean 1/p
inline int randgeo(double p) { return int(.999999 + log(randfrac())/log(1-p)); }

// Filter database stuff

extern int ADDRLEN; // set in db_generator.cc
#define ADDRBYTES ADDRLEN/8
#define MAXFILTERS 130000
#define MAXSTR 100
// #define NULL 0

struct filter {
  uint128_t sa; // IP source address
  uint128_t da; // IP destination address
  int sa_len; // IP source address mask length
  int da_len; // IP destination address mask length
  int sp[2]; // Transport source port range [low,high]
  int dp[2]; // Transport destination port range [low,high]
  int prot_num; // IP protocol field
  unsigned flags; // 16-bit flags field
  unsigned flags_mask; // 16-bit mask for flags
  int num_ext_field; // Number of extra header fields
  int *ext_field; // Pointer to array of extra header fields
};

// Do a deep copy of orig filter
inline void copy_filter(struct filter& copy, struct filter orig) {
  copy.sa            = orig.sa;
  copy.da            = orig.da;
  copy.sa_len        = orig.sa_len;
  copy.da_len        = orig.da_len;
  copy.sp[0]         = orig.sp[0];
  copy.sp[1]         = orig.sp[1];
  copy.dp[0]         = orig.dp[0];
  copy.dp[1]         = orig.dp[1];
  copy.prot_num      = orig.prot_num;
  copy.flags         = orig.flags;
  copy.flags_mask    = orig.flags_mask;
  copy.num_ext_field = orig.num_ext_field;
  if (copy.num_ext_field > 0) {
    copy.ext_field   = new int[copy.num_ext_field];
    for (int i = 0; i < copy.num_ext_field; i++) {
      copy.ext_field[i] = orig.ext_field[i];
    }
  } else {
    copy.ext_field   = NULL;
  }
}

struct range {
  int low;
  int high;
};

struct ppair {
  int slen;
  int dlen;
};

inline double factorial(int x) {
  double factx = 1;
  for (int i = 1; i <= x; i++) factx *= i;
  return factx;
}
