; ─────────────────────────── MAIN DEFAULT PARAMETERS ──────────────────────────
(define-param res     32)  ; resolution
(define-param kvecs   (list (vector3 0.234 0.45 0.321) (vector3 0.234 -0.45 0.321) )) ; list of k-vectors
(define-param rvecs   (list (vector3 -0.5 0.5 0.4) (vector3 0.5 -0.5 0.4) (vector3 0.5 0.5 -0.4) ))
(define-param nbands  4)   ; # of bands to compute
(define-param has-tr  #t)  ; whether assume B = 0 (#t) or apply apply a B-field along the cartesian z-direction (#f)
(define-param rad 0.25)
(define-param new-mesh-size 3)

; ──────────────────────────── PREP-WORK ────────────────────────────
; set calculation-specific parameters
(set! resolution res)
(set! num-bands  nbands)
(set! mesh-size  new-mesh-size)

; remove output of *.h5 unitcell files
(set! output-mu      (lambda () (print "... skipping output-mu\n")))
(set! output-epsilon (lambda () (print "... skipping output-epsilon\n")))

; ──────────────────── CRYSTAL SYSTEM/DIRECT BASIS ─────────────────────
(set! geometry-lattice (make lattice
    (size 1 1 1)
    (basis1 (list-ref rvecs 0))  (basis2 (list-ref rvecs 1))  (basis3 (list-ref rvecs 2))
    (basis-size (vector3-norm (list-ref rvecs 0)) (vector3-norm (list-ref rvecs 1)) (vector3-norm (list-ref rvecs 2)))
))

; ──────────────────── MATERIAL ─────────────────────
(define epsin 13)
(define epsin-diag    (vector3 14.00142849854971 14.00142849854971 13.0)) ; equiv. to a normalized B field of 0.4, applied to epsin=13
(define epsin-offdiag (cvector3 0.0+5.2i 0.0+0.0i 0.0+0.0i))

(define material-inside
    (cond
        (has-tr (make dielectric (epsilon epsin)) )
        (else   (make medium-anisotropic (epsilon-diag epsin-diag) (epsilon-offdiag epsin-offdiag)) )
    )
)
(define material-outside (make dielectric (epsilon 1)))

(use-modules (srfi srfi-1))  ; some additional fold-utilities; don't really need them, but nice to have
(define rad2 (* rad rad))
(define (f r)
    (let ((rc (lattice->cartesian r)))
        (- (abs (vector3-cdot rc rc)) rad2)
    )
)

(define (level-set-material r)
    (cond 
        ((< (f r) 0.0) material-inside)
        (else material-outside)
    )
)

(set! default-material (make material-function (material-func level-set-material)))

; ──────────────────── K-POINTS ─────────────────────
(set! k-points (map cartesian->reciprocal kvecs))

; ───────────────────────────────── DISPERSION ─────────────────────────────────
; run a dispersion calculation
(run-parity NO-PARITY true) ; equivalent to (run-all)

