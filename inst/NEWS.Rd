\name{NEWS}
\title{News for \R Package \pkg{wbacon}}
\encoding{UTF-8}
\section{CHANGES in wbacon VERSION 0.5-3 (2021-11-01)}{
    \subsection{Tweaks}{
        \itemize{
            \item In previous versions, the \code{summary} method for
                objects of class \code{wbaconlm} was implemented as a wrapper
                function of \code{stats::summary.lm}. But \code{summary.lm}
                detected that it has been called for an object other than a
                genuine \code{lm} object and issued a warning. This warning
                confused some users. Now we implemented a separate summary
                method.
        }
    }
}
\section{CHANGES in wbacon VERSION 0.5-2 (2021-10-07)}{
    \subsection{Misc}{
        \itemize{
            \item In calls to BLAS and LAPACK, we now use the same new
                macros that \code{R} uses in its code (\code{FCONE},
                \code{FC_LEN_T}, and \code{USE_FC_LEN_T}) since version 3.6.2.
                This approach is recommended in Writing R Extensions, see
                Chapter 6.6.1.
        }
    }
}
\section{CHANGES in wbacon VERSION 0.5-1 (2021-06-16)}{
    \subsection{BUG FIXES}{
        \itemize{
            \item fixed variadic macro ##__VA_ARGS__ because it is not
                portable (clang complains)
            \item fixed implicit type conversion of function 'med3' (only
                detected by clang on Apple M1), now int return type
            \item OpenMP issue: 'omp_set_num_threads' is now guarded by a
                conditional compilation directive
            \item OpenMP issue: removed '#pragma omp simd' directive in
                wbacon.c (line 480) because the optimizer was unable to
                perform the vectorization
        }
    }
}
\section{CHANGES in wbacon VERSION 0.5 (2021-06-13)}{
    \itemize{
        \item Initial release
    }
}
