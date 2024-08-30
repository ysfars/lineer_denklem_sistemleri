#lang racket
(require math/matrix)
(require math/vector)
(require math/array)

(define (vector-swap! vec i j)
  (let ((x (vector-ref vec i)))
    (vector-set! vec i (vector-ref vec j))
    (vector-set! vec j x)))

(define (vector-ref* x i j)(vector-ref (vector-ref x i) j))

(define (vector-map* katsayi satir) (vector-map (lambda (x) (* katsayi x)) satir))

;Yukarıdaki üç fonksiyon gauss-jordan ve cholesky fonksiyonlarında kullanılıyor.
          
(define (sms A) (square-matrix-size A)) ;; Kısaltma

(define (sonuc_yazdir A size)
  (for ([i size])
    (printf "x~a = ~a\n" i (matrix-ref A i 0))))
    
(define (det_2x2 A) (- (* (matrix-ref A 0 0) (matrix-ref A 1 1)) 
                       (* (matrix-ref A 0 1) (matrix-ref A 1 0))))

(define (kofaktor A) ;--> Determinant fonksiyonu.
  (if
   (equal? (sms A) 2) 
   (det_2x2 A)
   (for/sum ([i (sms A)])  
     (*(matrix-ref A 0 i) 
       (expt -1 (+ i 2))   
       (kofaktor
        (submatrix A       
                   (cdr(build-list (sms A) values))                 
                   (remove i (build-list (sms A) values))))))))

(define (cramer A b)
  (if (equal? (kofaktor A) 0)
      (error "Çözüm yok veya sonsuz çözüm var")
      (sonuc_yazdir
       (for*/matrix (sms A) 1 ([i (sms A)] [j 1])
         (/ (kofaktor
             (for*/matrix (sms A) (sms A) ([k (sms A)] [l (sms A)])
               (if (equal? i l)
                   (matrix-ref b k 0)
                   (matrix-ref A k l))))
            (kofaktor A)))
       (sms A))))

(define (gauss-jordan A b)
  (define rows (matrix->vector* A))
  (define B (matrix->vector* b))
  (for ([i (sms A)])
    (cond [(zero? (vector-ref* rows i i))
           (for ([j (in-range (add1 i) (sms A))]) ;--> Pivotlama
             #:final (not (zero? (vector-ref* rows j i)))
             (cond [(not (zero? (vector-ref* rows j i))) 
                    (vector-swap! rows i j)
                    (vector-swap! B i j)]))])

    (cond [(zero? (vector-ref* rows i i))
           (error "Çözüm yok veya sonsuz çözüm var")]) ;--> Determinant 0
    
    (vector-set! (vector-ref B i) 0
                 (/ (vector-ref* B i 0) (vector-ref* rows i i)))
    (vector-set! rows i (vector-map*                  ;--> Pivot elemanları kendisine bölünür.
                         (/ 1 (vector-ref* rows i i))
                         (vector-ref rows i))) 
    
    (for ([k (sms A)]) ;--> Diğer elemanların değeri 0 yapılır.
      (cond [(not (or (equal? i k) (equal? (vector-ref* rows k i) 0)))
             (vector-set! (vector-ref B k) 0
                          (- (vector-ref* B k 0) (* (vector-ref* rows k i) (vector-ref* B i 0))))
             (vector-set! rows k
                          (vector-map -
                                      (vector-ref rows k)
                                      (vector-map* (vector-ref* rows k i) (vector-ref rows i))))])))
  (sonuc_yazdir (vector*->matrix B) (sms A)))

;Pivotlama yapılmadığı için cholesky fonksiyonu sıfıra bölme hatası verebilir.
(define (cholesky A b)          
  (define U (matrix->vector* A))
  (define L (matrix->vector* (identity-matrix (sms A))))
  (define B (matrix->vector* b))
  (for ([i (sms A)]) 
    (for ([k (in-range (add1 i) (sms A))])
      (cond [(equal? (vector-ref* U k i) 0) (vector-set! (vector-ref L k) i 0)]
            [else
             (vector-set! (vector-ref L k) i (/ (vector-ref* U k i) (vector-ref* U i i)))
             (vector-set! U k
                          (vector-map
                           -
                           (vector-ref U k)
                           (vector-map* (/ (vector-ref* U k i) (vector-ref* U i i))
                                        (vector-ref U i))))])))
  ;Ax = LUx = LD = b --> Ux = D
  (for ([i (sms A)]) ;LD = b --> D bulunur ve B matrisine kopyalanır.
    (vector-set! (vector-ref B i) 0 (- (vector-ref* B i 0)
                                       (for/sum ([k i])
                                         (* (vector-ref* L i k) (vector-ref* B k 0))))))
  
  (for ([i (in-range (sub1 (sms A)) -1 -1)]) ;Ux = D --> x bulunur ve B Matrisine kopyalanır.
    (vector-set! (vector-ref B i) 0 (/ (- (vector-ref* B i 0)
                                          (for/sum ([k (in-range (add1 i) (sms A))])
                                            (* (vector-ref* U i k) (vector-ref* B k 0))))
                                       (vector-ref* U i i))))
  (printf "L matrisi: ~a\n\nU matrisi: ~a\n\n" L U)
  (sonuc_yazdir (vector*->matrix B) (sms A)))