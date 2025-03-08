//
//  DFrFT2D.swift
//  Cusp
//
//  Created by Daniel Clelland on 09/03/2025.
//

import Foundation
import Plinth

extension Matrix where Scalar == Double {
    
    public func dfrft2D(order: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).dfrft2D(order: order, setup: setup)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func dfrft2D(order: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        var result = self
        
        for i in 0..<shape.rows {
            result[row: i] = result[row: i].dfrft1D(order: order, setup: setup)
        }
        
        for j in 0..<shape.columns {
            result[column: j] = result[column: j].asRow().dfrft1D(order: order, setup: setup).asColumn()
        }
        
        return result
    }
    
}
