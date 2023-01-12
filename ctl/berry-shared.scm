; ----------------------------------------------------------------------
; ----- MISC SUB-ROUTINES SHARED AMONG BERRY PHASE IMPLEMENTATIONS -----
; ----------------------------------------------------------------------


; ----- MISC SUB-ROUTINES -----

; (check-integer-vec G tol) : checks whether a vector 'G' (vector3 struct.) equals (up to tolerance 'tol') its nearest 
;                             integer multiple vector; if not, throw an error specific to this calculation
(define (check-integer-vec G tol)
    (if (not (vector3-close? G (vector3 (round (vector3-x G)) (round (vector3-y G)) (round (vector3-z G))) tol))
        (begin
            (let ((error-msg (string-append "the loop considered in the (berry-phase ...) procedure "
                "must traverse an integer multiple of G vectors; currently, the loop's G vector is "
                "#(" (number->string (vector3-x G)) " " (number->string (vector3-y G)) " "
                (number->string (vector3-z G)) ") - i.e. not an integer multiple")))
                (error error-msg)
            )
        )
    )
)


; ----- MISC PRINT ROUTINES -----

; (print-berry-phases k-list phases) : print the print the berry phases and 
;                                      loop start k-point + G-vector along 
;                                      the loop (in reciprocal basis!)
(define (print-berry-phases k-list phases)
    (print "\nberry phases:")
    (for-each (lambda (j) (print ", " j)) (vector3->list (first k-list))) ; first k-point in loop (coords along G1, G2, G3)
    (for-each (lambda (j) (print ", " j)) (vector3->list (vector3- (last k-list) (first k-list)))); (integer) "vector of the loop" along G1, G2, G3
    (for-each (lambda (x) (for-each (lambda (y) (print ", " y)) x)) phases) ; print berry phases of bands/multiplets in concatenated style
    (print "\n\n")
)

; (print-berry-phase-output-struct bands-list) : print the output ordering to help the user 
;                                                understand the output/input relation
(define (print-berry-phase-output-struct bands-list)
    (print "\n\n---- Calculation of Berry phases by Wilson loop scheme ----\n")
    (print     "----       rows: loops; columns: bands/multiplets      ----\n\n")
    (print "berry phases:, loop-start k1, loop-start k2, loop-start k3, loop G1, loop G2, loop G3")
    (for-each (lambda (x) (cond ((= (length x) 1) (print ", band " x))
                                ((> (length x) 1) (print ", multiplet " x))  )
              ) bands-list)
    (print "\n\n")
)
