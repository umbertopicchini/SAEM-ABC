// XSum.c
// XSUM - Sum with error compensation
// The accuracy of the sum of floating point numbers is limited by the
// truncation error. E.g. SUM([1e16, 1, -1e16]) replies 0 instead of 1. The
// error grows with the length of the input, e.g. the error of SUM(RANDN(N, 1))
// is about EPS*(N / 10).
// Kahan, Knuth, Dekker, Ogita and Rump (and others) have derived some
// methods to reduce the influence of rounding errors, which are implemented
// here as fast C-Mex.
//
// Y = XSum(X, N, Method)
// INPUT:
//   X: Double array of any size.
//   N: Dimension to operate on. Optional, default: first non-singelton
//      dimension. Use the empty matrix [] to use this default together with
//      a 3rd argument.
//   Method: String: 'Double', 'Long', 'Kahan', 'Knuth', 'KnuthLong', 'Knuth2'.
//      Optional, default: 'Knuth'.
//      Not case sensitive. A description of the methods follows.
//   Call XSum without inputs to show the compilation date and availability of
//   long double methods.
//
// OUTPUT:
//   Y: Double array, equivalent to SUM, but with compensated error depending
//      on the Method. The high-precision result is rounded to double precision.
//      Y has the same size as X except for the N'th dimension, which is 1.
//      If X is empty, the empty matrix [] is replied.
//
// METHODS:
// The speed is compared to a single-threaded SUM for 1E6 elements. SUM is much
// faster with multi-threading and less than 1E4 elements, so only the relations
// between the methods have an absolute meaning.
// The accuracy depends on the input, so the stated values are estimations and
// must be seen as rules om thumb!
// Double: A thread-safe implementation of Matlab's SUM. At least in Matlab
//         2008a to 2009b the results of the multi-threaded SUM can differ
//         slightly from call to call.
//         Speed:    0.5 of time for SUM (MSVC++ 2008).
//         Accuracy: If the compilers supports accumulation in a 80 bit register
//                   (e.g. Open Watcom 1.8), 3 additional digits are gained.
//                   Otherwise the result is equivalent to SUM.
// Long:   The sum is accumulated in a 80 bit long double, if the compiler
//         supports this (e.g. LCC v3.8, Intel compilers).
//         Speed:    40% slower than SUM.
//         Accuracy: 3-4 more valid digits compared to SUM.
// Kahan:  The local error is subtracted from the next element. See [4].
//         Speed:    10% slower than SUM.
//         Accuracy: 1 to 3 more valid digits than SUM.
// Knuth:  Calculate the sum as if it is accumulated in a 128 bit float. This is
//         suitable for the most real world problems. See [1], [2], [3].
//         Speed:    As fast as SUM (MSVC++ 2008).
//         Accuracy: About 15 more valid digits than SUM.
// Knuth2: As Knuth, but the error is accumulated with a further level of
//         compensation. This is equivalent to an accumulation in a 196 bit
//         float. See [2] and INTLAB.
//         Speed:    60% slower than SUM (MSVC++ 2008).
//         Accuracy: 30 more valid digits than SUM.
// KnuthLong: As Knuth, but the temporary variables are long doubles, if this
//         is supported by the compiler.
//         Speed:    250% slower than SUM (LCC 3.8).
//         Accuracy: 21 more valid digits than SUM.
//
// NOTES: I've built a version using 384 bit QFLOATs with the free LCC v3.8
//   compiler (not the v2.4 shipped with Matlab): 104 digits but only 2%
//   speed of SUM. Feel free to mail me if you need instructions.
//
//   The subfunction DoubleDim1/N/V are a nice example of how to implement
//   functions operating on subvectors in arbitrary dimensions. This is much
//   faster and memory-efficient than using PERMUTE as in GRADIENT.
//
//   If anybody needs it, I publish XSum with multi-thread support, which
//   operates on complete subvectors separately. The idea of calculating
//   partial sums of large vectors in different threads is numerically unstable,
//   because the same input can yield different results!
//
//   Sorting the input increases the accuracy *sometimes*, but this is slow:
//   sorting 1E7 elements needs 30 times longer than adding them. For e.g.
//   SUM(RANDN(1, 1E7)) sorting reduces the accuracy! So this is no alternative.
//
// COMPILE:
//   mex -O XSum.c
// Linux: Consider C99 comments (thanks Sebastiaan Breedveld):
//   mex -O CFLAGS="\$CFLAGS -std=C99" XSum.c
// Pre-compiled MEX files (LCC v3.8, with Long Double, without QFloat):
//   http://www.n-simon.de/mex
//
// REFERENCES:
// [1] D.E. Knuth: "The Art of Computer Programming: Seminumerical
//   Algorithms", volume 2. Addison Wesley, Reading, Massachusetts,
//   second edition, 1981.
//
// [2] Takeshi Ogita and Siegfried M. Rump and Shin'ichi Oishi:
//   "Accurate Sum and Dot Product with Applications"
//   Proceedings of 2004 IEEE International Symposium on
//   Computer Aided Control Systems Design, Taipei, pages 152155, 2004
//   NOTE: Error in PDF "OgRuOi04a.pdf", TwoSum algorithm:
//     y = fl((a-(x-z))+(b + z))  ==>  y = fl((a-(x-z))+(b - z))
//
// [3] T.J. Dekker: "A Floating-Point Technique for Extending the
//   Available Precision". Numerische Mathematik, 18:224-242, 1971.
//
// [4] Linnainmaa, S.: Analysis of Some Known Methods of Improving the Accuracy
//   of Floating-point Sums, BIT 14, (1974), 167-202.
//
// Tested: Matlab 6.5, 7.7, 7.8, WinXP, [UnitTest]
//         Compilers: BCC5.5, LCC2.4/3.8, Open Watcom 1.8, MSVC++ 2008.
// Author: Jan Simon, Heidelberg, (C) 2009-2010 matlab.THISYEAR(a)nMINUSsimon.de
//
// See also: SUM.
// FEX: Alain Barraud: SUMBENCH (#23198), FORWARD_STABLE_SOLVER (#10668),
//      SUMDOTPACK (#8765)
// INTLAB: Prof. Dr. Siegfried M. Rump, http://www.ti3.tu-harburg.de/rump/intlab

/*
% $JRev: R0m V:012 Sum:m22wu8hRQASt Date:28-Feb-2010 02:39:51 $
% $License: BSD (use, copy, distribute on own risk, mention the author) $
% $File: Published\XSum\XSum.c $
% History:
% 001: 12-Feb-2010 21:51, First version: Sum with error compensations.
%      I built a version with QFloats of LCC v3.8.
%      XDot is under construction...
*/

#include "mex.h"
#include <string.h>
#include <stdlib.h>
#include <float.h>

// Assume 32 bit array dimensions for Matlab 6.5:
// See MEX option "compatibleArrayDims" for MEX in Matlab >= 7.7.
#ifndef MWSIZE_MAX
#define mwSize  int32_T           // Defined in tmwtypes.h
#define mwIndex int32_T
#define MWSIZE_MAX MAX_int32_T
#endif

// Prototypes:
void hello(void);
mwSize GetStep(const mwSize *XDim, const mwSize N);
mwSize FirstNonSingeltonDim(const mwSize *Xdim, const mwSize Xndim);

void KnuthDim1(double *X, const mwSize nX, const mwSize rX, double *Y);
void KnuthDimN(double *X, const mwSize *Xdim, const mwSize nX,
               const mwSize N, double *Y);
double KnuthV(double *X, double *XEnd, mwSize Step);

void Knuth2Dim1(double *X, const mwSize nX, const mwSize rX, double *Y);
void Knuth2DimN(double *X, const mwSize *Xdim, const mwSize nX,
                const mwSize N, double *Y);
double Knuth2V(double *X, double *XEnd, mwSize Step);

void KnuthLDim1(double *X, const mwSize nX, const mwSize rX, double *Y);
void KnuthLDimN(double *X, const mwSize *Xdim, const mwSize nX,
                const mwSize N, double *Y);
double KnuthLV(double *X, double *XEnd, mwSize Step);

void KahanDim1(double *X, const mwSize nX, const mwSize rX, double *Y);
void KahanDimN(double *X, const mwSize *Xdim, const mwSize nX,
               const mwSize N, double *Y);
double KahanV(double *X, double *XEnd, mwSize Step);

void DoubleDim1(double *X, const mwSize nX, const mwSize rX, double *Y);
void DoubleDimN(double *X, const mwSize *Xdim, const mwSize nX,
                const mwSize N, double *Y);

void LongDim1(double *X, const mwSize nX, const mwSize rX, double *Y);
void LongDimN(double *X, const mwSize *Xdim, const mwSize nX,
              const mwSize N, double *Y);

void SetPrecision(int Command);

// Different methods:
typedef enum {Unknown_SUM, Double_SUM, Long_SUM, Kahan_SUM,
              Knuth_SUM, KnuthLong_SUM, Knuth2_SUM} Method_t;

// Length of the method names (longest + 1):
#define Method_LEN 10

// Flag to show the warning for missing LONG support once only:
static int NoLong_ShowWarning = 1;

// Main function ===============================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mwSize *Xdim;
  mwSize       nX, *Ydim, Xndim, N;
  int          opOnDim1;
  Method_t     Method = Knuth_SUM;
  char         Method_In[Method_LEN];
  
  // Show general information:
  if (nrhs == 0) {
     hello();
     return;
  }
  
  // Check number of inputs and outputs:
  if (nrhs > 3) {
     mexErrMsgTxt("*** XSum[mex]: 1 to 3 inputs required.");
  }
  if (nlhs > 1) {
     mexErrMsgTxt("*** XSum[mex]: 1 output allowed.");
  }
  
  // Input must be a double array - is it worth to implement a sum with higher
  // accuracy for SINGLEs ?!
  if (!mxIsDouble(prhs[0])) {
     mexErrMsgTxt("*** XSum[mex]: X must be a double array.");
  }
  
  // Get first argument:
  nX    = mxGetNumberOfElements(prhs[0]);
  Xdim  = mxGetDimensions(prhs[0]);
  Xndim = mxGetNumberOfDimensions(prhs[0]);
  if (nX == 0) {   // Return the empty matrix if the input is empty:
     plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
     return;
  }
  
  // Get 2nd input [Dim] - default is the first non-singelton dimension:
  if (nrhs >= 2) {
     switch (mxGetNumberOfElements(prhs[1])) {
        case 0:  // Input [Dim] is empty - first non-singelton dimension:
           N        = FirstNonSingeltonDim(Xdim, Xndim);
           opOnDim1 = 1;
           break;
           
        case 1:
           // mxGetScalar needs a numeric input:
           if (!mxIsNumeric(prhs[1])) {
              mexErrMsgTxt("*** XSum[mex]: 2nd input [Dim] must be a "
                           "scalar double or empty.");
           }
           
           N = ((mwSize) mxGetScalar(prhs[1])) - 1;  // Zero based index
           if (N >= Xndim || N < 0) {
              mexErrMsgTxt("*** XSum[mex]: Dimension N out of range.");
           }
           
           // The fast algorithm for the 1st dimension works also for [1 x M]
           // row vectors. [1 x 1 x M] arrays are not caught.
           opOnDim1 = (N == 0 || (Xndim == 2 && N == 1 && *Xdim == 1));
           break;
           
        default:
           mexErrMsgTxt("*** XSum[mex]: 2nd input [Dim] must be "
                        "a scalar double or empty.");
     }
     
  } else {  // nrhs == 1: Operate on 1st non-singelton dimension:
     N        = FirstNonSingeltonDim(Xdim, Xndim);
     opOnDim1 = 1;
  }
  
  // The sum over a scalar is a scalar:
  if (Xdim[N] == 1) {
     plhs[0] = mxDuplicateArray(prhs[0]);
     return;
  }
  
  // Create output:
  if ((Ydim = (mwSize *)mxMalloc(Xndim * sizeof(mwSize))) == NULL) {
     mexErrMsgTxt("*** XSum[mex]: Cannot get memory for ndim.");
  }
  memcpy(Ydim, Xdim, Xndim * sizeof(mwSize));
  Ydim[N] = 1;
  plhs[0] = mxCreateNumericArray(Xndim, Ydim, mxDOUBLE_CLASS, mxREAL);
  mxFree(Ydim);
  
  // Parse 3rd input, default is KNUTH:
  if (nrhs == 3) {
     if (!mxIsChar(prhs[2])) {
        mexErrMsgTxt("*** XSum[mex]: 3rd input must be a string.");
     }
     
     mxGetString(prhs[2], Method_In, Method_LEN);
     if (strnicmp(Method_In, "Knuth", Method_LEN) == 0) {
        Method = Knuth_SUM;   // The default
     } else if (strnicmp(Method_In, "Knuth2", Method_LEN) == 0) {
        Method = Knuth2_SUM;
     } else if (strnicmp(Method_In, "KnuthLong", Method_LEN) == 0) {
        if (sizeof(long double) > sizeof(double)) {
           Method = KnuthLong_SUM;
        } else {              // Warn once if long doubles are doubles:
           if (NoLong_ShowWarning == 1) {
              mexWarnMsgIdAndTxt("JSim:XSum:NoLongDouble",
                "Using Knuth2 instead of Knuth_Long: "
                "XSum compiled without long double support!");
              NoLong_ShowWarning = 0;
           }
           Method = Knuth2_SUM;
        }
     } else if (strnicmp(Method_In, "Double", Method_LEN) == 0) {
        Method = Double_SUM;
     } else if (strnicmp(Method_In, "Long", Method_LEN) == 0) {
        Method = Long_SUM;
     } else if (strnicmp(Method_In, "Kahan", Method_LEN) == 0) {
        Method = Kahan_SUM;
     } else {             // Bad name, caught later
        Method = Unknown_SUM;
     }
  }
  
  // Call the core calculator for the specific method:
  switch (Method) {
     case Knuth_SUM:
        if (opOnDim1) {
           KnuthDim1(mxGetPr(prhs[0]), nX, Xdim[N], mxGetPr(plhs[0]));
        } else {
           KnuthDimN(mxGetPr(prhs[0]), Xdim, nX, N, mxGetPr(plhs[0]));
        }
        break;
        
     case Knuth2_SUM:
        if (opOnDim1) {
           Knuth2Dim1(mxGetPr(prhs[0]), nX, Xdim[N], mxGetPr(plhs[0]));
        } else {
           Knuth2DimN(mxGetPr(prhs[0]), Xdim, nX, N, mxGetPr(plhs[0]));
        }
        break;
        
     case KnuthLong_SUM:
        if (opOnDim1) {
           KnuthLDim1(mxGetPr(prhs[0]), nX, Xdim[N], mxGetPr(plhs[0]));
        } else {
           KnuthLDimN(mxGetPr(prhs[0]), Xdim, nX, N, mxGetPr(plhs[0]));
        }
        break;

     case Kahan_SUM:
        if (opOnDim1) {
           KahanDim1(mxGetPr(prhs[0]), nX, Xdim[N], mxGetPr(plhs[0]));
        } else {
           KahanDimN(mxGetPr(prhs[0]), Xdim, nX, N, mxGetPr(plhs[0]));
        }
        break;
        
     case Double_SUM:
        if (opOnDim1) {
           DoubleDim1(mxGetPr(prhs[0]), nX, Xdim[N], mxGetPr(plhs[0]));
        } else {
           DoubleDimN(mxGetPr(prhs[0]), Xdim, nX, N, mxGetPr(plhs[0]));
        }
        break;

     case Long_SUM:
        if (opOnDim1) {
           LongDim1(mxGetPr(prhs[0]), nX, Xdim[N], mxGetPr(plhs[0]));
        } else {
           LongDimN(mxGetPr(prhs[0]), Xdim, nX, N, mxGetPr(plhs[0]));
        }
        break;
        
     default:
        mexErrMsgTxt("*** XSum[mex]: 3rd input [Method] not recognized.");
  }
  
  return;
}

// =============================================================================
mwSize FirstNonSingeltonDim(const mwSize *Xdim, const mwSize Xndim)
{
  // Get first non-singelton dimension - zero based.
  mwSize N;
  
  for (N = 0; Xdim[N] == 1 && N < Xndim; N++) ;  // empty loop
  
  // No non-singelton dim found - reply 0:
  if (N == Xndim) {
     N = 0;
  }
  
  return (N);
}

// =============================================================================
mwSize GetStep(const mwSize *Xdim, const mwSize N)
{
  // Get step size between elements of a subvector in the N'th dimension.
  const mwSize *XdimEnd, *XdimP;
  mwSize       Step;
  
  Step    = 1;
  XdimEnd = Xdim + N;
  for (XdimP = Xdim; XdimP < XdimEnd; Step *= *XdimP++) ; // empty loop
  
  return (Step);
}

// =============================================================================
void KnuthDim1(double *X, const mwSize nX, const mwSize rX, double *Y)
{
  // Calculate Knuth's sum over vectors in the 1st dimension.
  
  // METHOD:
  //   % Sum(X): See reference [2]
  //     n = length(X);
  //     for i = 2:n
  //       [X(i), X(i-1)] = TwoSum(X(i), X(i-1))
  //     end
  //     Sum = sum(X(1:n-1)) + X(n);
  //
  //   % 6 ops, no branch better to optimize, [1]:
  //   function1 [x, y] = TwoSum(a, b)
  //     x = fl(a + b)
  //     z = fl(x - a)
  //     y = fl((a - (x - z)) + (b - z))
  //
  //   % 4 ops, branching impedes pipelining in the processor, [3]:
  //   function2 [x, y] = TwoSum(a, b)
  //     x = fl(a + b)
  //     if abs(a) >= abs(b)
  //       y = fl(b - (x - a))
  //     else
  //       y = fl(a - (x - b))
  
  double *XR, *XEnd, h, q, x, z;
  
  XEnd = X + nX;
  while (X < XEnd) {
     XR = X + rX;
     h  = *X++;
     q  = 0.0;
     while (X < XR) {
        x  = h + *X;                        // Predictor
        z  = x - h;                         // Corrector
        q += ((h - (x - z)) + (*X++ - z));  // Accumulate the local error
        h  = x;
     }
     *Y++ = h + q;                          // Add accumulated local errors
  }
  
  return;
}

// -----------------------------------------------------------------------------
void KnuthDimN(double *X, const mwSize *Xdim, const mwSize Xlen,
               const mwSize N, double *Y)
{
  // Calculate Knuth's sum over vectors in N'th dimension.
  double *XEnd, *SliceEnd;
  mwSize Step, SumStep;
  
  // Get step size between elements of the subvector and its extent:
  Step    = GetStep(Xdim, N);
  SumStep = Step * Xdim[N];
  
  // Loop over all dimension:
  XEnd = X + Xlen;
  while (X < XEnd) {
     for (SliceEnd = X + Step; X < SliceEnd; X++) {
        // X goes through the elements of the subvector in the N'th dimension:
        *Y++ = KnuthV(X, X + SumStep, Step);
     }  // end for
     
     // Move pointer to the next slice:
     X += SumStep - Step;
  }
 
  return;
}

// -----------------------------------------------------------------------------
double KnuthV(double *X, double *XEnd, mwSize Step)
{
  // Calculate Knuth's sum for a subvector.
  // [Step] is the distance between elements of the subvector.
  double q, h, x, z;

  h  = *X;
  X += Step;
  q  = 0.0;
  while (X < XEnd) {
     x  = h + *X;
     z  = x - h;
     q += ((h - (x - z)) + (*X - z));
     h  = x;
     X += Step;
  }

  return (h + q);
}

// =============================================================================
void Knuth2Dim1(double *X, const mwSize nX, const mwSize rX, double *Y)
{
  // Calculate Knuth's sum over vectors in the 1st dimension, accumulate the
  // local errors with one extra stage of correction. Finally the sum is
  // obtained as if it was accumulated in a 196 bit double.
  
  double *XR, *XEnd, h, q1, q2, x, z, x2, z2, e1;
  
  XEnd = X + nX;
  while (X < XEnd) {
     XR = X + rX;
     h  = *X++;
     q1  = 0.0;
     q2  = 0.0;
     while (X < XR) {
        x = h + *X;
        z = x - h;
        
        // Accumulate the local error with a further level of correction:
        e1  = ((h - (x - z)) + (*X++ - z));
        x2  = q1 + e1;
        z2  = x2 - q1;
        q2 += ((q1 - (x2 - z2)) + (e1 - z2));
        q1  = x2;
        
        h = x;
     }
     *Y++ = h + q1 + q2;
  }
  
  return;
}

// -----------------------------------------------------------------------------
void Knuth2DimN(double *X, const mwSize *Xdim, const mwSize Xlen,
                const mwSize N, double *Y)
{
  // Calculate Knuth's sum over vectors in N'th dimension.
  double *XEnd, *SliceEnd;
  mwSize Step, SumStep;
  
  // Get step size between elements of the subvector and its extent:
  Step    = GetStep(Xdim, N);
  SumStep = Step * Xdim[N];
  
  // Loop over all dimension:
  XEnd = X + Xlen;
  while (X < XEnd) {
     for (SliceEnd = X + Step; X < SliceEnd; X++) {
        // X goes through the elements of the subvector in the N'th dimension:
        *Y++ = Knuth2V(X, X + SumStep, Step);
     }  // end for
     
     // Move pointer to the next slice:
     X += SumStep - Step;
  }
 
  return;
}

// -----------------------------------------------------------------------------
double Knuth2V(double *X, double *XEnd, mwSize Step)
{
  // Calculate Knuth's sum over vectors in the 1st dimension, accumulate the
  // local errors with one extra stage of correction. Finally the sum is
  // obtained as if it was accumulated in a 196 bit double.
  // [Step] is the distance between elements of the subvector.
  double h, q1, q2, x, z, x2, z2, e1;
  
  h   = *X;
  X  += Step;
  q1  = 0.0;
  q2  = 0.0;
  while (X < XEnd) {
     x = h + *X;
     z = x - h;
     
     // Accumulate the local error with a further level of correction:
     e1  = ((h - (x - z)) + (*X - z));
     x2  = q1 + e1;
     z2  = x2 - q1;
     q2 += ((q1 - (x2 - z2)) + (e1 - z2));
     q1  = x2;
     
     h  = x;
     X += Step;
  }
  
  return (h + (q1 + q2));
}

// =============================================================================
void KnuthLDim1(double *X, const mwSize nX, const mwSize rX, double *Y)
{
  // Calculate Knuth's sum over vectors in the 1st dimension with intermediate
  // values stored in long doubles.
  // If the compiler uses doubles for long doubles, the error correction fails
  // completely due to setting the precision to 64 bit mantissa. Then the
  // result are comparable with the uncompensated SUM!
  
  double *XR, *XEnd;
  long double h, q, x, z;
  
  SetPrecision(64);
  
  XEnd = X + nX;
  while (X < XEnd) {
     XR = X + rX;
     h  = (long double) *X++;
     q  = 0.0L;
     while (X < XR) {
        x  = h + *X;
        z  = x - h;
        q += ((h - (x - z)) + (*X++ - z));
        h  = x;
     }
     *Y++ = (double) (h + q);
  }
  
  SetPrecision(0);
  return;
}

// -----------------------------------------------------------------------------
void KnuthLDimN(double *X, const mwSize *Xdim, const mwSize Xlen,
                const mwSize N, double *Y)
{
  // Calculate Knuth's sum over vectors in N'th dimension with long doubles.
  // This is identical to KnuthDimN, just the subfunction for processing the
  // subvectors differs.
  
  double *XEnd, *SliceEnd;
  mwSize Step, SumStep;
  
  SetPrecision(64);
  
  // Get step size between elements of the subvector and its extent:
  Step    = GetStep(Xdim, N);
  SumStep = Step * Xdim[N];
  
  // Loop over all dimension:
  XEnd = X + Xlen;
  while (X < XEnd) {
     for (SliceEnd = X + Step; X < SliceEnd; X++) {
        // X goes through the elements of the subvector in the N'th dimension:
        *Y++ = KnuthLV(X, X + SumStep, Step);
     }  // end for
     
     // Move pointer to the next slice:
     X += SumStep - Step;
  }
 
  SetPrecision(0);
  return;
}

// -----------------------------------------------------------------------------
double KnuthLV(double *X, double *XEnd, mwSize Step)
{
  // Calculate Knuth's sum for a subvector with long doubles.
  // [Step] is the distance between elements of the subvector.
  // Caller set precision to 64 bits mantissa.
  
  long double q, h, x, z;
  
  h  = (long double) *X;
  X += Step;
  q  = 0.0L;
  while (X < XEnd) {
     x  = h + *X;
     z  = x - h;
     q += (h - (x - z)) + (*X - z);
     h  = x;
     X += Step;
  }
  
  return ((double) (h + q));
}

// =============================================================================
void KahanDim1(double *X, const mwSize nX, const mwSize rX, double *Y)
{
  // Calculate Kahan's sum over vectors of the 1st dimension.
  
  double *XR, *XEnd, c, t, x, s;
  
  // Kahan sum does not work with 64 bit precision, because this conceals the
  // cancellation error!!! 53 bit precision is the default in Matlab, but the
  // user could have changed this by "system_dependent('setprecision', 64)".
  
  // Method used in XBLAS:DASUM (see www.netlib.org):
  //   % This is a double compenstated sum
  //   high = 0.0;
  //   low  = 0.0;
  //   for i = 1:length(X) {
  //      X_elem = X(i);
  //      t1   = high + X_elem;
  //      e    = t1 - high;
  //      t2   = ((X_elem - e) + (high - (t1 - e))) + low;
  //      high = t1 + t2;
  //      low  = t2 - (high - t1);
  //    }
  //    Sum = high;
  
  SetPrecision(53);
  
  XEnd = X + nX;
  while (X < XEnd) {
     XR = X + rX;
     s  = *X++;
     c  = 0.0;
     while (X < XR) {
       x = *X++ - c;
       t = s + x;
       c = (t - s) - x;
       s = t;
     }
     
     *Y++ = s - c;
  }
  
  SetPrecision(0);
  return;
}

// -----------------------------------------------------------------------------
void KahanDimN(double *X, const mwSize *Xdim, const mwSize Xlen,
               const mwSize N, double *Y)
{
  // Calculate Kahan's sum over vectors in N'th dimension.
  double *XEnd, *SliceEnd;
  mwSize Step, SumStep;
  
  SetPrecision(53);  // Not working with 64!
    
  // Get step size between elements of the subvector and its extent:
  Step    = GetStep(Xdim, N);
  SumStep = Step * Xdim[N];

  // Loop over all dimension:
  XEnd = X + Xlen;
  while (X < XEnd) {
     for (SliceEnd = X + Step; X < SliceEnd; X++) {
        // X goes through the elements of the subvector in the N'th dimension:
        *Y++ = KahanV(X, X + SumStep, Step);
     }  // end for
     
     // Move pointer to the next slice:
     X += SumStep - Step;
  }
 
  SetPrecision(0);
  return;
}

// -----------------------------------------------------------------------------
double KahanV(double *X, double *XEnd, mwSize Step)
{
  // Calculate Kahan's sum for a subvector.
  // [Step] is the distance between elements of the subvector.
  double c = 0.0, s, t, x;
  
  // Cheap first step:
  s  = *X;
  X += Step;
  
  // Loop over elements 2 to end:
  while(X < XEnd) {
    x  = *X - c;
    t  = s + x;
    c  = (t - s) - x;
    s  = t;
    X += Step;
  }
  
  return (s);
}

// =============================================================================
void DoubleDim1(double *X, const mwSize nX, const mwSize rX, double *Y)
{
  // Calculate sum over vectors of the 1st dimension in double precision.
  
  double *XR, *XEnd;
  double Sum;
  
  SetPrecision(64);

  XEnd = X + nX;
  while (X < XEnd) {
     Sum = 0.0;
     for (XR = X + rX; X < XR; Sum += *X++) ;  // empty loop
     *Y++ = Sum;
  }
  
  SetPrecision(0);
  return;
}

// -----------------------------------------------------------------------------
void DoubleDimN(double *X, const mwSize *Xdim, const mwSize Xlen,
                const mwSize N, double *Y)
{
  // Calculate sum over vectors in N'th dimension with long double precision.
  double *XEnd, *SliceEnd, *V, *VEnd;
  mwSize Step, SumStep;
  double Sum;
  
  SetPrecision(64);
  
  // Get step size between elements of the subvector and its extent:
  Step    = GetStep(Xdim, N);
  SumStep = Step * Xdim[N];
  
  // Loop over all dimension:
  XEnd = X + Xlen;
  while (X < XEnd) {
     for (SliceEnd = X + Step; X < SliceEnd; X++) {
        // V goes through the elements of the subvector in the N'th dimension:
        Sum  = 0.0L;
        VEnd = X + SumStep;
        for (V = X; V < VEnd; V += Step) {
           Sum += *V;
        }
        *Y++ = (double) Sum;
     }  // end for
     
     // Move pointer to the next slice:
     X += SumStep - Step;
  }
 
  SetPrecision(0);
  return;
}

// =============================================================================
void LongDim1(double *X, const mwSize nX, const mwSize rX, double *Y)
{
  // Calculate sum over vectors of the 1st dimension in long double precision.
  //
  // This algorithm can benefit from the 64 bit precision even if the compiler
  // uses doubles as long doubles: for the accumulation [Sum] is stored in a
  // 80 bit register. This happens e.g. for Open Watcom 1.8, but not for LCC
  // v2.4, which is shipped with Matlab, or MSVC++ 2008.
  
  double *XR, *XEnd;
  long double Sum;
  
  SetPrecision(64);

  XEnd = X + nX;
  while (X < XEnd) {
     Sum = 0.0L;
     for (XR = X + rX; X < XR; Sum += *X++) ;  // empty loop
     *Y++ = (double) Sum;
  }
  
  SetPrecision(0);
  return;
}

// -----------------------------------------------------------------------------
void LongDimN(double *X, const mwSize *Xdim, const mwSize Xlen,
              const mwSize N, double *Y)
{
  // Calculate sum over vectors in N'th dimension with long double precision.
  double *XEnd, *SliceEnd, *V, *VEnd;
  mwSize Step, SumStep;
  long double Sum;
  
  SetPrecision(64);
  
  // Get step size between elements of the subvector and its extent:
  Step    = GetStep(Xdim, N);
  SumStep = Step * Xdim[N];
  
  // Loop over all dimension:
  XEnd = X + Xlen;
  while (X < XEnd) {
     for (SliceEnd = X + Step; X < SliceEnd; X++) {
        // V goes through the elements of the subvector in the N'th dimension:
        Sum  = 0.0L;
        VEnd = X + SumStep;
        for (V = X; V < VEnd; V += Step) {
           Sum += *V;
        }
        *Y++ = (double) Sum;
     }  // end for
     
     // Move pointer to the next slice:
     X += SumStep - Step;
  }
 
  SetPrecision(0);
  return;
}

// =============================================================================
void SetPrecision(int Command)
{
  // Set floating point environment for Intel processor.
  // The long double calculations work with 80 bit numbers only, if the
  // precision is set to 64 bits. With 53 bits, the calculations have double
  // precision only.
  // The Kahan sum fails when processed with 64 bit precision, because this
  // steals the necessary bits from the correction of the local error.
  //
  // It seems that Matlab resets the floating point processor flags at the
  // standard exit of a Mex function. But at least Matlab 5.3 and 6.5 leave the
  // the flags untouched on exits through mexErrMsgTxt.
  //
  // TESTED FOR BCC AND LCC ONLY!
  // I have no idea how this works on processors, which are not supported by
  // the function _control87 !
  //
  // INPUT:
  //   Command: int, precision of floating poiunt calculations.
  //            Allowed values: 24 (single), 53 (double), 64 (long double),
  //                            0 (reset to previous value)
  
#if defined(__LCC__)           // Leading underscores needed for macros:
#define MCW_PC _MCW_PC
#define PC_24  _PC_24
#define PC_53  _PC_53
#define PC_64  _PC_64
#endif

  static unsigned int PC_Back = 0;
  
  if (Command == 0) {                   // Restore original settings
     _control87(PC_Back, MCW_PC);
  } else {
     PC_Back = _control87(MCW_PC, 0);   // Backup original settings
     if (Command == 24) {
        _control87(PC_24, MCW_PC);      // Single precision
     } else if (Command == 53) {
        _control87(PC_53, MCW_PC);      // Double precision
     } else if (Command == 64) {
        _control87(PC_64, MCW_PC);      // Full precision
     } else {
        mexErrMsgTxt("*** XSum[mex]: Bad command in SetPrecision.");
     }
  }
  
  return;
}

// =============================================================================
void hello(void)
{
  // Show some information.
  // Most important: availability of LONG DOUBLE
  
  char *hasLong, *yes = "yes", *no = "no";
  
  if (sizeof(long double) > sizeof(double)) {
     hasLong = yes;
  } else {
     hasLong = no;
  }
  
  // Actually the compiler would be important also, but it is not available
  // as macro.
  mexPrintf("XSum: Sum with compensated errors\n"
            "Author: Jan Simon\n"
#if defined(__DATE__) && defined(__TIME__)
            "Compiled: " __DATE__ " " __TIME__ "\n"
#endif
            "Long double: %s\n", hasLong);
  
  return;
}
