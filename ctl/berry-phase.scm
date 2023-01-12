; DESCRIPTION: Defines a function (berry-phase k-list bands-list) that calculates 
;              Abelian/non-Abelian (i.e. single/multiband) Berry phases via Wilson-loops.
; INPUT: k-list, a list of k-points (a closed, contractible or noncontractible, loop)
;                Eg.: (interpolate 20 (list (vector3 kx -0.5 0) (vector3 kx 0.5 0)))
;                     a loop with fixed kx value (=kx), varying along ky in 20 
;                     equidistant steps
;        bands-list, list of list of bands. Band-indices in inner lists need not be 
;                incremening or consecutive/neighboring (though this is more well-tested)
;                Eg.: (list (list 1) (list 2 3) (list 4) (list 5 6 8)) computes the 
;                     Abelian Berry phases of bands 1 and 4, and the non-Abelian Berry 
;                     phases of the (2,3) and (5,6,8) multiplets.
; OUTPUT: list of lists of (possibly non-Abelian) Berry phases, in the order of bands-list 
;         for loops along k-list.
;
; AUTHORS: This implementation was written by Thomas Christensen (Jan, 2018 -- June, 2019)
;          The inspiration for this implementation was Ling Lu's much simpler and more
;          minimalistic implementation (Jan 12, 2017; https://github.com/LingLuIOP/mpb).
;          The generalizations made here reflects a desire to: 
;             1. Avoid repeated k-solve calls in the typical case where several multiplets
;                are requested, 
;             2. To allow non-consecutive (and even non-incrementing) choices of bands in 
;                a multiplet,
;             3. Facilitate basic printing of output and input parameters

(include "aux.scm")
(include "berry-shared.scm")

(define (berry-phase k-list bands-list)
    (set! num-bands (maximum (flatten bands-list))) ; set max number of bands based on bands-list
    ; temporary variables defined in letstatement:
    ; - ev,     previous eigenvector saved
    ; - uk,     list of real-space wavefunctions at (H field) first k-point
    ; - ukG,    list of real-space wavefunctions at (B field) last k-point
    ; - W,      product of sqmatrices
    ; - G-loop, the difference in k after a closed loop
    ; - consec, list of booleans indicating whether the bands in bands-list are consecutive/neighboring
    ; - phases, actual berry phases (output)
    (let ( (ev (list)) (uk (list)) (ukG (list)) (W (list)) (G-loop (vector3)) (consec (list)) (phases (list)) )
        ; get G-vec associated with the loop
        (set! G-loop (vector3- (last k-list) (first k-list) ) ) 
        (check-integer-vec G-loop 0.0001) ; make sure G-loop is an integer multiple of G vectors 

        ; store the list of eigenfields of the first k-point
        (solve-kpoint (first k-list)) ; solve for first element of k-list
        (map (lambda (bands) 
                (let ((ukt (list)))
                    (map (lambda (band)
                            (get-hfield band) ; take the hfield (not the bfield; distinction matters when case mu \noteq 1)
                            (set! ukt (cons (field-copy cur-field) ukt))) 
                        bands)
                    (set! uk (cons (reverse ukt) uk))
                ) 
                (set! ev (cons (get-eigenvectors (minimum bands) (+ (- (maximum bands) (minimum bands)) 1)) ev )) ; eigenvectors in list (may contain more, if non-consecutive)
                (set! W (cons (sqmatrix-diagm (ones-list (length bands))) W)) ; initialized as unit-matrices (list of sqmatrix'es)
                (set! consec (append consec (list (consecutive? bands)))) ; whether the set of bands is consecutive or not (list of booleans)
             )  bands-list)
        (set! ev (reverse ev))  (set! W (reverse W))  (set! uk (reverse uk)) ; iterated cons gives reverse ordering - undo that here

        
        ; INNER-PRODUCT OF STATES ON THE K-LIST W = <k0|k1><k1|k2>...<k(N-1)|kN><kN|k0*e(-iGr)>           [the e(-iGr) part is done
        ;                                         = <Bk0|Hk1><Bk1|Hk2>...<Bk(N-1)|HkN><BkN|Hk0*e(-iGr)>    separately below]
        (do ((i 1 (+ i 1))) ((= i (length k-list))) ; loop over remaining k-points
            (solve-kpoint (list-ref k-list i))
            (do ((bi 0 (+ bi 1))) ((= bi (length bands-list))) ; loop over bands-list ("bands indices": bi)
                (let ((bands (list-ref bands-list bi)))
                    ; dot product does dot(ev, state_j) = <B(ev)|H(state_j)>, i.e. the overlap
                    (let* ((overlap (dot-eigenvectors (list-ref ev bi) ; bands from previous k-point
                                                      (minimum bands))) ; bands from current k-point
                           (sub-overlap (cond ((list-ref consec bi) overlap)
                                              (else (pick-sub-parts overlap (zero-base bands)) ) 
                                        ))
                          )  
                          (list-set! W bi (sqmatrix-mult (list-ref W bi) sub-overlap)) 
                    )
                    (list-set! ev bi ; update bands to be used as previous k-point in next k-point
                               (get-eigenvectors (minimum bands)                             ; lowest band
                                                 (+ (- (maximum bands) (minimum bands)) 1) ) ; band-span
                    )
                )
            )
        )
 
        ; store the list of eigenfields at the last k-point (shouldn't be necessary, but ... oh well)
        (map (lambda (bands) 
                (let ((ukt (list)))
                    (map (lambda (band)
                            (get-bfield band) ; take the bfield (not the hfield; distinction matters when case mu \noteq 1)
                            (set! ukt (cons (field-copy cur-field) ukt))) 
                        bands)
                    (set! ukG (cons (reverse ukt) ukG))
                ) 
             ) bands-list)
        (set! ukG (reverse ukG)) ; again, iterated cons needs reordering

        ; FIRST/LAST (=EQUAL) EIGENVECTORS CONDITIONS
        ; the boundary conditions between the initial and final k points is enforced via the mutual G-phased overlap 
        (do ((bi 0 (+ bi 1))) ((= bi (length bands-list))) ; loop over bands-list 
            ; compute Gphase: the phase differences between first & last k-point
            (let ((Gphase (list))) ; empty Gphase list for each new set of bands/multiplets
                (do ((bii 0 (+ bii 1))) ((= bii (length (list-ref bands-list bi)) )) ; loop over lists in bands-list
                    (set! Gphase (cons (integrate-fields (lambda (r u1 u2) 
                                                                (*  (exp (* (vector3-dot G-loop r) (* 2 pi 0-1i))) 
                                                                    (vector3-cdot u1 u2))  ) ; conjugated dot-prod: u1* dot u2
                                                         (list-ref (list-ref ukG bi) bii)   ; field at last k-point
                                                         (list-ref (list-ref uk  bi) bii) ) ; field at first k-point
                                       Gphase) )
                )   
                 
                ; BUILD A DIAGONAL SQMATRIX HERE, THEN MULTIPLY W, THEN RETURN THE PHASES OF EIGENVALS
                (list-set! W bi (sqmatrix-mult (list-ref W bi) (sqmatrix-diagm (reverse Gphase)))) ; enforce the boundary cond. on W[bi] and reorder Gphase cf. iterated cons
            )
            ; berry phases: obtained as phases of eigvals of W[bi] matrix (sorted, ascending) 
            (set! phases (cons (sort (map (lambda (x) (- (imag-part (log x))))
                                          (sqmatrix-eigvals (list-ref W bi)))  <)
                                phases) ) ; stored as list of lists  
        )
        (set! phases (reverse phases)) ; reorder cf. iterated cons

        (print-berry-phases k-list phases) ; print results
        (begin phases) ; force output (return value of set! is unspecified)
    )
)

; ----- SQMATRIX SUB-ROUTINE -----

; (pick-sub-parts A lst) : returns an sqmatrix that is a submatrix of the 
;                          sqmatrix 'A', including only the rows and columns 
;                          corresponding to those in the (incrementing) list
;                          'lst' of zero-based indices.
;      Example: (pick-sub-parts (sqmatrix-diagm (list 1 2 3 4)) (list 1 3))
;               returns the matrix [2 0; 0 4] 
(define (pick-sub-parts A lst)
    (let* ((n (length lst)) (sub-A (sqmatrix-diagm (zeros-list n))))
        (map (lambda (i)
                (map (lambda (j)
                        (sqmatrix-set sub-A i j (sqmatrix-ref A (list-ref lst i) (list-ref lst j))))
                    (iota n)))
            (iota n))
        (begin sub-A)
    )
)
