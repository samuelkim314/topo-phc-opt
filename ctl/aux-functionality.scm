
; (vector3->list x) : convert a vector3 object 'r' to an ordinary list; this is 
;                     occasionally needed for (print ...), (for-each ...), or (map ...)
(define (vector3->list r)
    (list (vector3-x r) (vector3-y r) (vector3-z r))
)

; (bool->int x)     : converts a boolean x (= #f or #t) to an integer (0 or 1). 
;                     occassionally useful for arithmetic methods on bools
(define (bool->int x)
    (case x
        ((#f) 0)
        ((#t) 1)
    )
)

; (map-bool->int X) : converts a list of booleans X to a list of integers
(define (map-bool->int X)
    (map (lambda (x) (bool->int x)) X)
)

; --- FIND DIMENSIONALITY OF LATTICE ---
; (dim-lattice)    : returns the dimensionality of the lattice, based on the 
;                    number of no-size (=1e-20) entries in 'size property of 
;                    geometry-lattice
(define (lattice-size)
    (vector3->list (object-property-value geometry-lattice 'size)) 
)
(define (lattice-nosize)
    (map (lambda (x) (<= x 1e-20)) (lattice-size))
)
(define (dim-lattice)
    (- 3 (fold-left + 0 (map-bool->int (lattice-nosize))))
)

; (time expr N) : timing utility function: prints the average time taken (in seconds) to execute
;                 the stump expression `expr` (of the type `(lambda () ...)`) `N` times.
;                 Example:      (time (lambda () (sleep 1)) 10) 
;                               Prints: 10.0034527 "s"
(define time (lambda (expr N)
    (let ((t1 (get-internal-real-time)))
        (do ((i 1 (+ i 1))) ((> i N))
            (expr)
        )
        (display (exact->inexact (/ (- (get-internal-real-time) t1) (* N internal-time-units-per-second))))
        (display " s\n")
    )
))


; (keep-from-idx) : keep all the elements of the list `l` from `idx` to its end (assuming 
;                   1-based indexing).
;                   Example: (keep-from-idx (list 1 2 3 4 5 6) 4) => (4 5 6)
(define (keep-from-idx l idx)
    (cond 
        ((<= idx 1) (begin l))
        (else       (keep-from-idx (cdr l) (- idx 1)))
    )
)
