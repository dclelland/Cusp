//
//  FrFT2D.swift
//  Cusp
//
//  Created by June Russell on 02/03/2025.
//

import Foundation
import Numerics
import Plinth

//extension Matrix where Scalar == Float {
//    
//    public func frft2D(order: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
//        return ComplexMatrix(real: self).frft2D(order: order, setup: setup)
//    }
//    
//}
//
//extension ComplexMatrix where Scalar == Float {
//    
//    public func frft2D(order: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
//        var result = self
//        
//        for i in 0..<shape.rows {
//            result[row: i] = result[row: i].frft1D(order: order, setup: setup)
//        }
//        
//        for j in 0..<shape.columns {
//            result[column: j] = result[column: j].frft1D(order: order, setup: setup)
//        }
//        
//        return result
//    }
//    
//}

extension Matrix where Scalar == Double {
    
    public func frft2D(order: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).frft2D(order: order, setup: setup)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func frft2D(order: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        var result = self
        
        for i in 0..<shape.rows {
            result[row: i] = result[row: i].frft1D(order: order, setup: setup)
        }
        
        for j in 0..<shape.columns {
            result[column: j] = result[column: j].frft1D(order: order, setup: setup)
        }
        
        return result
    }
    
}
