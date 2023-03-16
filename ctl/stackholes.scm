; input values
(define-param rvecs    (list (vector3 1 0 0) (vector3 0 1 0) (vector3 0 0 1)))
(define-param mat1ff 0.5) ; filling fraction of the 1st material
(define-param hole-r 0.1)
(define-param epsin1    12.0)
(define-param epsin2    4.0)
(define-param epsout    1.0)


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

; ────────────────────── SET MATERIAL PROPERTIES ───────────────────────
(set! default-material (make dielectric (epsilon epsin1)))
(define mat2 (make dielectric (epsilon epsin2)))
(set! geometry
    (list
        (make block (size (vector3-norm (list-ref rvecs 0)) 
                        (vector3-norm (list-ref rvecs 1))
                        (* (- 1 mat1ff) (vector3-norm (list-ref rvecs 2)))))
        (make cylinder (radius hole-r) (height (vector3-norm (list-ref rvecs 2))))
    ))

