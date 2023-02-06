; ----- AUXILIARY FUNCTIONS -----

(define pi (* 4 (atan 1))) ; 3.14159...

; (maximum lst) : finds largest element of a (flat) list 'lst'
;               by just "splatting" 'lst' into the arguments of max
(define (maximum lst)
  (apply max lst)
)

; (minimum lst) : finds largest element of a (flat) list 'lst'
;               by just "splatting" 'lst' into the arguments of max
(define (minimum lst)
  (apply min lst)
)

; (zero-base lst) : given a list 'lst', returns a corresponding lst whose elements are
;                   shifted relative to the minimum element of 'lst', such that the
;                   minimum element of 'lst' is sent to 0. Useful for relative indexing.
;      Example: (zero-base (list 4 7 9)) returns (0 3 5)
(define (zero-base lst)
    (map (lambda (x) (- x (minimum lst))) lst) ; "zeroed" version of 'lst'; minimum element runs from 0
)

; (vector3->list x) : convert a vector3 object 'r' to an ordinary list; this is 
;                     occasionally needed for (print ...), (for-each ...), or (map ...)
(define (vector3->list r)
    (list (vector3-x r) (vector3-y r) (vector3-z r))
)

; (flatten lst) : flattens a list of lists or list of pairs or mixed
;                 list of pairs and lists 'lst' to a single "flat list"
;                 recursive implementation [https://stackoverflow.com/a/8387641]
(define (flatten lst)
    (cond ((null? lst) '())
        ((pair? lst) (append (flatten (car lst)) (flatten (cdr lst))))
        (else (list lst)))
)

; (last lst) : last element of a list 'lst' [the "opposite" of (first lst)]
;              can also be implemented as: 
;                   (define (last lst) 
;                       (first (list-tail lst (- (length lst) 1))) 
;                   )
(define (last lst)
    (cond   ((null? (cdr lst)) (first lst)) ; if only 1 element in list, return this element
            (else (last (cdr lst)))
    )
)

; (list-set lst j val) : sets the value 'val' at the index 'j' in a list 'lst'
;                        conceptually, this is equivalent to Julia/Matlab/Python'ish syntax
;                            lst[j] .= val
;     [note that lists are zero-indexed in Scheme; https://stackoverflow.com/a/7382392]
(define (list-set! lst j val)
    (if (zero? j)
        (set-car! lst val)
        (list-set! (cdr lst) (- j 1) val))
)

; (append-val-until lst val n) : appends 'val' to a list 'lst' 'n' times
;                                e.g. (append-val-until (list 1 2) 4 2) gives (list 1 2 4 4)
(define (append-val-until lst val n)
    (cond ((> n 0) (append-val-until (cons val lst) val (- n 1)))
          (else lst))
)

; (ones-list n) : create a list of 1 elements of length n
;                 e.g. (ones-list 5) is ( 1 1 1 1 1 )
(define (ones-list n) (append-val-until (list) 1 n))

; (zeros-list n) : create a list of 0 elements of length n
;                  e.g. (zeros-list 5) is ( 0 0 0 0 0 )
(define (zeros-list n) (append-val-until (list) 0 n))

; (consecutive? lst) : returns true if a list 'lst' contains only consecutive, incrementing elements, e.g.:
;                      (1 2 3) returns true; (1 3 4) or (2 1 3) returns false. 
;                      A single-element is defined as consecutive.
(define (consecutive? lst)
    (cond ((null? (cdr lst)) #t)
          (else (cond ((equal? (+ (car lst) 1) (car (cdr lst))) (consecutive? (cdr lst)))
                      (else #f)       )))
)

; (sqmatrix-print A) : prints an sqmatrix object, taking no special care to align columns
(define (sqmatrix-print A)
    (let ((size (sqmatrix-size A)))
        (map (lambda (i)
            (map (lambda (j)
                (print (sqmatrix-ref A i j) ", ")
            ) (iota size))
            (print "\n")
        ) (iota size))
    )
)