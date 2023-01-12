; modules 
(use-modules (srfi srfi-1))   ; additional fold/reduce-utilities
(use-modules (srfi srfi-4))   ; for f64-vector functionality
(use-modules (ice-9 vlist))   ; for vlists

; includes
(include "../aux-functionality.scm")
; (define (vector3 x y z) (vector x y z)
; (define (vector3* a v) (vector (* a (vector-ref v 0)) (* a (vector-ref v 1)) (* a (vector-ref v 2)) ))

; variables
(define uc-gvecs (list (vector3 0 0 0) (vector3 0 -1 -1) (vector3 0 -1 1) (vector3 -1 0 -1) (vector3 -1 0 1) (vector3 1 0 -1) (vector3 1 0 1) (vector3 0 1 -1) (vector3 0 1 1) (vector3 -1 -1 0) (vector3 1 1 0) (vector3 1 -1 0) (vector3 -1 1 0) (vector3 0 0 -2) (vector3 0 0 2) (vector3 0 -1 -2) (vector3 0 -1 2) (vector3 -1 0 -2) (vector3 -1 0 2) (vector3 1 0 -2) (vector3 1 0 2) (vector3 0 1 -2) (vector3 0 1 2) (vector3 0 -2 0) (vector3 -2 0 0) (vector3 2 0 0) (vector3 0 2 0) (vector3 -1 -1 -2) (vector3 -1 -1 2) (vector3 1 1 -2) (vector3 1 1 2) (vector3 1 -1 -2) (vector3 1 -1 2) (vector3 -1 1 -2) (vector3 -1 1 2) (vector3 0 -2 -1) (vector3 0 -2 1) (vector3 -2 0 -1) (vector3 -2 0 1) (vector3 2 0 -1) (vector3 2 0 1) (vector3 0 2 -1) (vector3 0 2 1) (vector3 0 -2 -2) (vector3 0 -2 2) (vector3 -2 0 -2) (vector3 -2 0 2) (vector3 2 0 -2) (vector3 2 0 2) (vector3 0 2 -2) (vector3 0 2 2)))
(define uc-coefs (list 0.006509988713824982-0.07430991596602736i -0.14858501604673846+0.037798932181882614i 0.14858501604673843-0.03779893218188265i 0.1485850160467384-0.037798932181882586i -0.1485850160467385+0.0377989321818826i -0.1485850160467385+0.03779893218188263i 0.1485850160467385-0.03779893218188265i 0.14858501604673843-0.03779893218188259i -0.1485850160467385+0.03779893218188263i -0.05796045325126349+0.009834853427071608i -0.057960453251263505+0.00983485342707161i -0.07397410695454079-0.040613880878154165i -0.07397410695454078-0.040613880878154186i -0.055001995038864486+0.07342833704372212i -0.055001995038864486+0.07342833704372212i -0.045988480619645444+0.00967774700948442i 0.045988480619645444-0.009677747009484424i -0.045988480619645444+0.009677747009484417i 0.04598848061964546-0.00967774700948442i 0.04598848061964546-0.009677747009484422i -0.04598848061964546+0.009677747009484426i 0.04598848061964546-0.009677747009484419i -0.04598848061964547+0.009677747009484422i 0.0006733988398186295+0.005833383133757143i 0.0006733988398186281+0.005833383133757143i 0.0006733988398186298+0.0058333831337571435i 0.0006733988398186282+0.0058333831337571435i 0.022141415260758895-0.020356404143648046i 0.02214141526075885-0.02035640414364806i 0.022141415260758898-0.02035640414364805i 0.022141415260758853-0.02035640414364806i -0.008413912285000654-0.014572480450529158i -0.008413912285000658-0.014572480450529165i -0.008413912285000656-0.014572480450529165i -0.008413912285000666-0.014572480450529179i 0.0005102563342371311-0.029259171073349383i 0.0005102563342371369-0.029259171073349393i -0.0005102563342371481+0.029259171073349383i -0.0005102563342371405+0.029259171073349383i -0.0005102563342371329+0.029259171073349383i -0.0005102563342371369+0.029259171073349393i 0.0005102563342371483-0.029259171073349393i 0.0005102563342371407-0.029259171073349393i 0.002588281223189957-0.0006278232264738109i 0.002588281223189956-0.0006278232264738123i 0.002588281223189957-0.0006278232264738104i 0.002588281223189958-0.0006278232264738094i 0.002588281223189958-0.000627823226473811i 0.002588281223189957-0.000627823226473811i 0.002588281223189958-0.0006278232264738104i 0.002588281223189957-0.0006278232264738113i))

(define twopi 6.283185307179586)
(define uc-coefs-real  (map (lambda (coef) (real-part coef))    uc-coefs)) ; hoist some operations outside r-loop
(define uc-coefs-imag  (map (lambda (coef) (imag-part coef))    uc-coefs))
(define uc-gvecs-twopi (map (lambda (g)    (vector3*  twopi g)) uc-gvecs))

; default testing point
(define rp    (vector3 .2 .1 .3))
(define f64rp (f64vector (vector-ref rp 0) (vector-ref rp 1) (vector-ref rp 2)))
; -------------------------------------------------------------

(define (vector3-dot-manual v w) 
    (+ (* (vector-ref v 0) (vector-ref w 0))
       (* (vector-ref v 1) (vector-ref w 1))
       (* (vector-ref v 2) (vector-ref w 2)))
)

(define (f64vector3-dot-manual v w) 
    (+ (* (f64vector-ref v 0) (f64vector-ref w 0))
       (* (f64vector-ref v 1) (f64vector-ref w 1))
       (* (f64vector-ref v 2) (f64vector-ref w 2)))
)

(define (f64vector3-dot-manual-mixed v w) 
    (+ (* (f64vector-ref v 0) (vector-ref w 0))
       (* (f64vector-ref v 1) (vector-ref w 1))
       (* (f64vector-ref v 2) (vector-ref w 2)))
)

(define dot-test             (lambda () (vector3-dot rp rp)))
(define dot-test-manual      (lambda () (vector3-dot-manual rp rp)))
(define f64f-test-dot-manual (lambda () (f64vector3-dot-manual f64rp f64rp)))
(time dot-test             100000)
(time dot-test-manual      100000)
(time f64f-test-dot-manual 100000)

; -------------------------------------------------------------

(define (calc-test-dot-manual r)
    (fold + 0.0
        (map (lambda (g) (vector3-dot-manual g r)
            ) uc-gvecs-twopi)
    )
)

(define (calc-test-dot r)
    (fold + 0.0
        (map (lambda (g) (vector3-dot g r)
            ) uc-gvecs-twopi)
    )
)

(define f-dot-manual (lambda () (calc-test-dot-manual rp)))
(define f-dot (lambda () (calc-test-dot rp)))

(time f-dot-manual 1000)
(time f-dot 1000)

; -------------------------------------------------------------
; -------------------- IMPLEMENTATIONS ------------------------
; -------------------------------------------------------------

; original
(define (calc-fourier-sum r)
    (fold + 0.0
        (map (lambda (g coef-real coef-imag)
                (let ((g-dot-r (vector3-dot g r)))
                    (- (* coef-real (cos g-dot-r))
                       (* coef-imag (sin g-dot-r)))  )
            ) uc-gvecs-twopi uc-coefs-real uc-coefs-imag)
    )
)

; use manual dot function instead of vector3-dot
(define (calc-fourier-sum-manual r)
    (fold + 0.0
        (map (lambda (g coef-real coef-imag)
                (let ((g-dot-r (vector3-dot-manual g r)))
                    (- (* coef-real (cos g-dot-r))
                       (* coef-imag (sin g-dot-r)))  )
            ) uc-gvecs-twopi uc-coefs-real uc-coefs-imag)
    )
)

; two iterators (lists) instead of 3
(define (calc-fourier-sum-manual-a r)
    (fold + 0.0
        (map (lambda (g coef)
                (let ((g-dot-r (vector3-dot-manual g r)))
                    (- (* (real-part coef) (cos g-dot-r))
                       (* (imag-part coef) (sin g-dot-r)))
                )
            ) uc-gvecs-twopi uc-coefs)
    )
)

; one iterator only: list of vectors (also use reduce; doesn't matter though)
(define iter-vecs 
    (map (lambda (g coef)
            (vector (vector3-x g) (vector3-y g) (vector3-z g) (real-part coef) (imag-part coef))
        ) uc-gvecs-twopi uc-coefs)
)
(define (calc-fourier-sum-manual-b r)
    (reduce + 0.0
        (map (lambda (g-and-coefs)
                (let ((g-dot-r (vector3-dot-manual g-and-coefs r)))
                    (- (* (vector-ref g-and-coefs 3) (cos g-dot-r))
                       (* (vector-ref g-and-coefs 4) (sin g-dot-r)))
                )
            ) iter-vecs)
    )
)

; vlist
(define iter-vecs-as-vlist (list->vlist iter-vecs))
(define (calc-fourier-sum-manual-c r)
    (vlist-fold + 0.0
        (vlist-map (lambda (g-and-coefs)
                (let ((g-dot-r (vector3-dot-manual g-and-coefs r)))
                    (- (* (vector-ref g-and-coefs 3) (cos g-dot-r))
                       (* (vector-ref g-and-coefs 4) (sin g-dot-r)))
                )
            ) iter-vecs-as-vlist)
    )
)

; vector of vectors
(define iter-vec-vecs (list->vector iter-vecs)) ; benefit: stored contiguously in memory, effectively like an array
(define N (vector-length iter-vec-vecs))
(define Nm1 (- N 1))
(define (calc-fourier-sum-manual-d r)
    (do ((i 0 (+ 1 i))
         (v 0.0 (+ v (let* ((g-and-coefs (vector-ref iter-vec-vecs i))
                                (g-dot-r (vector3-dot-manual g-and-coefs r)))
                                (- (* (vector-ref g-and-coefs 3) (cos g-dot-r))
                                   (* (vector-ref g-and-coefs 4) (sin g-dot-r)))                
                     )))
        )
        ((>= i N) v)
    )
)

(define (calc-fourier-sum-manual-e r)
    (let ((x 0.0))
        (do ((i 0 (+ 1 i)))
            ((>= i N))
                (let* ( (g-and-coefs (vector-ref iter-vec-vecs i))
                        (g-dot-r (vector3-dot-manual g-and-coefs r)))
                        (set! x (+ x (- (* (vector-ref g-and-coefs 3) (cos g-dot-r))
                                     (* (vector-ref g-and-coefs 4) (sin g-dot-r)))  ))        
                )
        )
        (begin x)
    )
)

; manual licm of the components of r 
; (versions f through h have identical performance; a matter of taste which is preferred)
(define (calc-fourier-sum-manual-f r)
    (let ((v 0.0) (rx (vector-ref r 0)) (ry (vector-ref r 1)) (rz (vector-ref r 2)))
        (do ((i 0 (+ 1 i)))
            ((>= i N))
                (let* ((g-and-coefs (vector-ref iter-vec-vecs i))
                       (g-dot-r     (+ (* (vector-ref g-and-coefs 0) rx)
                                       (* (vector-ref g-and-coefs 1) ry)
                                       (* (vector-ref g-and-coefs 2) rz)))
                       )
                        (set! v (+ v (- (* (vector-ref g-and-coefs 3) (cos g-dot-r))
                                        (* (vector-ref g-and-coefs 4) (sin g-dot-r)))  ))        
                )
        )
        (begin v)
    )
)

(define (calc-fourier-sum-manual-g r)
    (let ( (rx (vector-ref r 0))   ; licm of referencing r components
           (ry (vector-ref r 1))
           (rz (vector-ref r 2)) )
        (do ((i 0 (+ 1 i))
             (v 0.0 (+ v (let* ( (g-and-coefs (vector-ref iter-vec-vecs i))
                                 (gx          (vector-ref g-and-coefs 0))
                                 (gy          (vector-ref g-and-coefs 1))
                                 (gz          (vector-ref g-and-coefs 2))
                                 (c-real      (vector-ref g-and-coefs 3))
                                 (c-imag      (vector-ref g-and-coefs 4))
                                 (g-dot-r     (+ (* gx rx) (* gy ry) (* gz rz))) )
                                ; add real part of exponential term to v
                                (- (* c-real (cos g-dot-r)) (* c-imag (sin g-dot-r)))             
                         )))
            )
            ((>= i N) v) ;  continue until test is true (then return v
        )
    )
)

(define (calc-fourier-sum-manual-h r)
    (let ( (v 0.0)               ; loop summand
           (rx (vector-ref r 0)) ; licm of referencing r components
           (ry (vector-ref r 1)) 
           (rz (vector-ref r 2)) ) 
        (do ((i 0 (+ i 1))) ; loop counter
            ((>= i N) v)    ; continue until test is true (then return v)
            ; loop body
            (let* ( (g-and-coefs (vector-ref iter-vec-vecs i))
                    (gx          (vector-ref g-and-coefs 0))
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

; f64
; (figured using f64 would give some gain: doesn't appear to)
(define iter-f64vecs
    (map (lambda (g coef)
            (f64vector (vector3-x g) (vector3-y g) (vector3-z g) (real-part coef) (imag-part coef))
        ) uc-gvecs-twopi uc-coefs)
)
(define iter-f64vecs-as-vlist (list->vlist iter-f64vecs))
(define (f64calc-fourier-sum-manual r)
    (vlist-fold + 0.0
        (vlist-map (lambda (g-and-coefs)
                (let ((g-dot-r (f64vector3-dot-manual-mixed g-and-coefs r)))
                    (- (* (f64vector-ref g-and-coefs 3) (cos g-dot-r))
                       (* (f64vector-ref g-and-coefs 4) (sin g-dot-r)))
                )
            ) iter-f64vecs-as-vlist)
    )
)

; ----------------------------------------------------------
; -------------------- BENCHMARKING ------------------------
; ----------------------------------------------------------

(define f-fourier-orig     (lambda () (calc-fourier-sum          rp)))
(define f-fourier-manual   (lambda () (calc-fourier-sum-manual   rp)))
(define f-fourier-manual-a (lambda () (calc-fourier-sum-manual-a rp)))
(define f-fourier-manual-b (lambda () (calc-fourier-sum-manual-b rp)))
(define f-fourier-manual-c (lambda () (calc-fourier-sum-manual-c rp)))
(define f-fourier-manual-f64 (lambda () (f64calc-fourier-sum-manual rp)))
(define f-fourier-manual-d (lambda () (calc-fourier-sum-manual-d rp)))
(define f-fourier-manual-e (lambda () (calc-fourier-sum-manual-e rp)))
(define f-fourier-manual-f (lambda () (calc-fourier-sum-manual-f rp)))
(define f-fourier-manual-g (lambda () (calc-fourier-sum-manual-g rp)))
(define f-fourier-manual-h (lambda () (calc-fourier-sum-manual-h rp)))

(define bench (lambda (n)
    (print "Original :\t")  (time f-fourier-orig       n)
    (print "Manual :\t")    (time f-fourier-manual     n)
    (print "Manual-a:\t")   (time f-fourier-manual-a   n)
    (print "Manual-b:\t")   (time f-fourier-manual-b   n)
    (print "Manual-c:\t")   (time f-fourier-manual-c   n)
    (print "Manual-f64:\t") (time f-fourier-manual-f64 n)
    (print "Manual-d:\t")   (time f-fourier-manual-d   n)
    (print "Manual-e:\t")   (time f-fourier-manual-e   n)
    (print "Manual-f:\t")   (time f-fourier-manual-f   n)
    (print "Manual-g:\t")   (time f-fourier-manual-g   n)
    (print "Manual-h:\t")   (time f-fourier-manual-h   n)
))

(define minibench (lambda (n)
    (print "Manual-f:\t")   (time f-fourier-manual-f   n)
    (print "Manual-g:\t")   (time f-fourier-manual-g   n)
    (print "Manual-h:\t")   (time f-fourier-manual-h   n)
))

(print "--- All benchmarks ---\n")
(bench 1000)

(print "\n--- Zoom'ed benchmarks ---\n")
(minibench 50000)

(quit)