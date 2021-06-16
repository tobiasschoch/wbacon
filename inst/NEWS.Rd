\name{NEWS}
\title{News for \R Package \pkg{wbacon}}
\encoding{UTF-8}
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
