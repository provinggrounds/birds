;;; Set-up a PA of the icosahedral lattice
;;
;; control the tolerance
(set! tolerance 1e-6)
(set! eigensolver-block-size -11)
(set-param! mesh-size 37)
(set-param! num-bands 350)
(define-param res 16) ; should be at least 8
(define-param res_x 16) ; should be a power of 2
(define-param res_y 16)

(define-param NN 299)
(define-param Chi "0.5")
(define-param Method "Bare")
(define-param LengthUCX (sqrt 299))
(define-param LengthUCY (sqrt 299))
(define-param LengthUCZ 1)

(set!  geometry-lattice  (make lattice 
			 (size  1 1 no-size) ; number of boxes
                         (basis-size LengthUCX LengthUCY LengthUCZ )
                         (basis1 1 0 0)
                         (basis2 0 1 0)
                         ))

;;; material parameters

(define pi (* 4 (atan 1)))

(set-param! epsSi 11.56) ; dielectric constant
(set-param! epsAir 1) ; dielectric constant

(define-param epsBack epsAir) ; dielectric constant
(define-param epsCylH epsSi) ; dielectric constant
(define-param epsCylV epsSi) ; dielectric constant

(define dielBack (make dielectric (epsilon epsBack)))
(define dielCylH (make dielectric (epsilon epsCylH)))
(define dielCylV (make dielectric (epsilon epsCylV)))


;;; Fill the space with background dielectric
(set! default-material  dielBack)

;; Set-up parameters that define the ball-and-stick PA icosahedral lattice

(define-param CylinderVerticalRadius   0.20)

(define-param theta 0)
(define-param phi   (* pi 0.5))

(define-param CV_R CylinderVerticalRadius)
;vertical cylinder heights and axes
(define CV_H infinity)
;(define CV_A (vector3 0 0 1)) ; A: axis to define orientation

(define CV_A (lambda (theta phi) (vector3 (* (sin theta) (cos phi) ) (* (sin theta) (sin phi) ) (* (cos theta) 1 ))))

(include "twospecies_define_r0.3000_r0.0500.ctl" )

(include "twospecies_make.ctl")

;;; Set-up the k-space analysis
;;Set-up the high-symmetry points of the Brillouin zone for the SC lattice,

(define Gamma     (vector3 0    0   0))     		; Gamma point
(define X         (vector3 0.5  0   0))	  	        ; X point
(define M         (vector3 0.5  0.5 0))			; M point
(define R         (vector3   0  0.5 0))        	        ; R point

;; Set up the k-space path for band structure: Gamma-X-M-R-Gamma
(define-param k-interp 3)
(set! k-points (interpolate k-interp (list  Gamma X M R Gamma))); region bigger than irred brillouin zone
;;; (set! k-points (list  Gamma X))

;;;; Calculate length of K space path segments
;; k-points are in rec. lattice units, convert in cartesian (there is a 2 \pi factor missing)
(define GammaC  (vector3* (* 2 3.14159)  (reciprocal->cartesian Gamma)))	; Gamma point
(define XC      (vector3* (* 2 3.14159)  (reciprocal->cartesian X)))           	; X point
(define MC      (vector3* (* 2 3.14159)  (reciprocal->cartesian M)))           	; M point
(define RC      (vector3* (* 2 3.14159)  (reciprocal->cartesian R)))     	; R point

(define GammaX 	(vector3-norm (vector3- GammaC XC)))
(define XM 	(vector3-norm (vector3- XC MC)))
(define MR 	(vector3-norm (vector3- MC RC)))
(define RGamma 	(vector3-norm (vector3- RC GammaC)))

;; Set-up the running parameters

(set-param! resolution (vector3 (* res_x res) (* res_y res) 1)) ; use a 8x8x8 grid

(define fname (string-append "_twospecies_MHUDS_r0.3000_r0.0500"))

;; Run calculation
;(run-tm);
(run-tm (output-at-kpoint X  fix-efield-phase output-efield-z)) ; run simulation and output E-field at point X
