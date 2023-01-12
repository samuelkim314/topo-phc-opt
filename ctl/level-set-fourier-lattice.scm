; input values
(define-param rvecs    (list (vector3 1 0 0) (vector3 0 1 0) (vector3 0 0 1)))
(define-param uc-gvecs (list (vector3 0 0 0))) ; initialized as constant function
(define-param uc-coefs (list 1.0+0.0i))
(define-param uc-level  0.0)
(define-param epsin    10.0)
(define-param epsout    1.0)
(define-param epsin-diag     '()) ; a vector3:  (εxx εyy εzz) [Cartesian axes]
(define-param epsin-offdiag  '()) ; a cvector3: (εxy εxz εyz) [Cartesian axes]
(define-param epsout-diag    '()) ; a vector3:  (εxx εyy εzz) [Cartesian axes]
(define-param epsout-offdiag '()) ; a cvector3: (εxy εxz εyz) [Cartesian axes]


; ──────────────────── CRYSTAL SYSTEM/DIRECT BASIS ─────────────────────
(cond 
    ((equal? dim 1) ; 1D
        (set! geometry-lattice (make lattice
            (size 1 no-size no-size) 
            (basis1 (list-ref rvecs 0))
            (basis-size (vector3-norm (list-ref rvecs 0)) 1 1)
        ))
    )
    ((equal? dim 2) ; 2D
        (set! geometry-lattice (make lattice
            (size 1 1 no-size) 
            (basis1 (list-ref rvecs 0)) (basis2 (list-ref rvecs 1))
            (basis-size (vector3-norm (list-ref rvecs 0)) (vector3-norm (list-ref rvecs 1)) 1)
        ))
    )
    ((equal? dim 3) ; 3D
        (set! geometry-lattice (make lattice
            (size 1 1 1)
            (basis1 (list-ref rvecs 0)) (basis2 (list-ref rvecs 1)) (basis3 (list-ref rvecs 2))
            (basis-size (vector3-norm (list-ref rvecs 0)) (vector3-norm (list-ref rvecs 1)) (vector3-norm (list-ref rvecs 2)))
        ))
    )
)

; ──────────────────────────── PERMITTIVITY ────────────────────────────
; calculate the real part of the fourier sum
;
; carefully optimized version of fourier summation: see `experiments/experiments-calc-fourier.scm` 
; for detailed benchmarking; we additionally play a manual `(load ...)` trick here to ensure the 
; method gets compiled. 
; the call to `(load ...)` defines a method `calc-fourier-sum-fast`, which is functionally equivalent
; to the out-commented `calc-fourier-sum` method below (i.e. returns the real part of the fourier sum)
(add-to-load-path ".") ; ensure the ctl directory is in load path (+ avoid issues w/ paths relative to current working dir)
(load "calc-fourier-sum-fast.scm")

; simple implementation (for reference & posterity): 
;       (use-modules (srfi srfi-1))  ; some additional fold-utilities; don't really need them, but nice to have
;       (define twopi 6.283185307179586)
;       (define uc-gvecs-twopi (map (lambda (g) (vector3* twopi g)) uc-gvecs)) ; hoist outside r-loop
;       (define (calc-fourier-sum r) ; I've spent a full day trying (many!) variants of this to improve its speed; 
;           (reduce + 0.0            ; but it doesn't seem doable when sticking to Scheme alone
;               (map (lambda (g coef)
;                       (let ((g-dot-r (vector3-dot g r)))
;                           (- (* (real-part coef) (cos g-dot-r))
;                              (* (imag-part coef) (sin g-dot-r)))  )
;                    ) uc-gvecs-twopi uc-coefs)
;           )
;       )

; define associated level-set function (< uc-level => inside; > uc-level => outside)
(define (level-set-eps r)
    (cond ((< (calc-fourier-sum-fast r) uc-level) epsin)
          (else epsout))
)

(define (level-set-tensor-eps r) ; only used if either `epsin-offdiag` or `epsout-offdiag` is not null; assumes equi-diagonal epsilon
    (cond ((< (calc-fourier-sum-fast r) uc-level)
              (make medium-anisotropic (epsilon-diag    epsin-diag)
                                       (epsilon-offdiag epsin-offdiag))   )
          (else
              (make medium-anisotropic (epsilon-diag    epsout-diag)
                                       (epsilon-offdiag epsout-offdiag))  )
    )
)

; ──────────────────────────── PERMEABILITY ────────────────────────────
; optionally, a non-uniform, non-unity permeability can be set [default: μ(r) = 1]
(define-param uc-gvecs-mu '())
(define-param uc-coefs-mu '())
(define-param uc-level-mu 0.0)
(define-param muin        1.0)
(define-param muout       1.0)

; same as above, but for permeability mu; defines calc-fourier-sum-mu from uc-gvecs-mu and uc-coefs-mu
(load "calc-fourier-sum-fast-mu-copy.scm")

; define associated level-set function (< uc-level => inside; > uc-level => outside)
(define (level-set-mu r)
    (cond ((< (calc-fourier-sum-fast-mu r) uc-level-mu) muin)
          (else muout))
)

; define joint epsilon+mu level set function
(define (level-set-eps-and-mu r)
    (make medium (epsilon (level-set-eps r)) (mu (level-set-mu r)) )
)

; ────────────────────── SET MATERIAL PROPERTIES ───────────────────────
(cond 
    ((null? uc-gvecs-mu)  ; nonuniform permittivity & uniform permeabilty (=1)
        (cond 
            ((and (null? epsin-offdiag) (null? epsout-offdiag))
                (set! default-material (make material-function (epsilon-func level-set-eps)))
            )
            (else         ; epsilon tensor
                (print "\nAnisotropic material: a tensorial permittivity with nonzero off-diagonal components was specified...\n")
                (set! default-material (make material-function (material-func level-set-tensor-eps)))
            )
        )
    )
    (else                 ; nonuniform permittivity _and_ permability
        (set! force-mu? true) ; # must set explicitly since we have no objects with μ ≠ 1 (see https://github.com/NanoComp/mpb/issues/39)
        (set! default-material (make material-function (material-func level-set-eps-and-mu)))
    )
)