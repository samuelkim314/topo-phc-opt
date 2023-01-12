(define (print-symeigs kidx W w eigs)
    (print "sym-eigs:")
    (print ", " kidx)
    (print ", (" (op->string W w) ")")
    (for-each (lambda (e) (print ", " e)) eigs)
    (print "\n")
)

(define (op->string W w)
    (let ( (Wlist (list (list (vector3-x (vector3-x W)) (vector3-x (vector3-y W)) (vector3-x (vector3-z W)))
                        (list (vector3-y (vector3-x W)) (vector3-y (vector3-y W)) (vector3-y (vector3-z W)))
                        (list (vector3-z (vector3-x W)) (vector3-z (vector3-y W)) (vector3-z (vector3-z W)))))
           (wlist (list (vector3-x w) (vector3-y w) (vector3-z w)))
           (dim (dim-lattice))
           (xyzt  "")
         )

        (do ((n 0 (+ n 1))) ((= n dim))      ; rows
            (do ((m 0 (+ m 1))) ((= m dim))  ; columns
                (let* ((Wnm (list-ref (list-ref Wlist n) m)) )
                    (if (not (zero? Wnm))
                        (set! xyzt (string-append xyzt
                                       (number->formatted-string Wnm)
                                       (if (= m 0) "x" (if (= m 1) "y" "z")) )
                        )
                    )
                )
            )
            (let ((wn (list-ref wlist n)))
                (if (not (zero? wn)) 
                    (set! xyzt (string-append xyzt 
                                              (if (positive? wn) "+" "")
                                              (number->string wn)))
                )
                (if (not (= n (- dim 1)))
                    (set! xyzt (string-append xyzt ","))
                )
            )
        )
        (begin xyzt) ; return xyzt
    )
)


(define (number->formatted-string x)
    (cond 
        ((zero? x) "")
        (else
            (string-append 
                (if (positive? x) "+" "-")
                (if (not (= (abs x) 1))
                        (number->string (if (integer? x) (inexact->exact (abs x)) (abs x)))
                        ""
                )
            )
        )
    )
)