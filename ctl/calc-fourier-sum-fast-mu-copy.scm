; This file is a _direct_ copy of `calc-fourier-sum-fast.scm` except that all variables
; and methods have `-mu` prepended to them. This is obviously awful, but it seemed the 
; easiest way to achieve the aims of `calc-fourier-sum-fast` without breaking the 
; solution to the tricky compilation issues, while still getting the equivalent for mu.

; Variables
(define (scm-vector3* a v) (vector (* a (vector-ref v 0)) (* a (vector-ref v 1)) (* a (vector-ref v 2)) )) ; defined solely so the file can be run from guile directly, instead of requiring mpb
(define uc-gvecs-mu-twopi (let ((twopi 6.283185307179586)) (map (lambda (g) (scm-vector3* twopi g)) uc-gvecs-mu)))

; Store uc-gvecs-mu-twopi and uc-coefs-mu contiguously in memory, as an N-vector of 5-vectors; 
; effectively, this is like an array
(define vec-of-vecs-g-and-coefs-mu
    (list->vector               ; N-vector
        (map (lambda (g coef)   ; 5-vector (gx, gy, gz, c-real, c-imag)
            (vector (vector-ref g 0) (vector-ref g 1) (vector-ref g 2) (real-part coef) (imag-part coef))
        ) uc-gvecs-mu-twopi uc-coefs-mu)
    )
)
(define uc-N-mu (vector-length vec-of-vecs-g-and-coefs-mu))

; Main function
(define (calc-fourier-sum-fast-mu r)
    (let ( (v 0.0)               ; loop summand
           (rx (vector-ref r 0)) ; licm of referencing r components
           (ry (vector-ref r 1)) 
           (rz (vector-ref r 2)) ) 
        ; no appropriate mapping routine for this vector-scenario: use do loop
        (do ((i 0 (+ i 1))) ; loop counter
            ((>= i uc-N-mu) v) ; continue until test is true (then return v)
            ; loop body
            (let* ( (g-and-coefs (vector-ref vec-of-vecs-g-and-coefs-mu i))
                    (gx          (vector-ref g-and-coefs 0)) ; w/ 2*pi factors
                    (gy          (vector-ref g-and-coefs 1))
                    (gz          (vector-ref g-and-coefs 2))
                    (c-real      (vector-ref g-and-coefs 3))
                    (c-imag      (vector-ref g-and-coefs 4))
                    (g-dot-r     (+ (* gx rx) (* gy ry) (* gz rz))) )
                ; add real part of exponential term to v
                (set! v (+ v (- (* c-real (cos g-dot-r)) 
                                (* c-imag (sin g-dot-r)) )))
            )
        )
    )
)
