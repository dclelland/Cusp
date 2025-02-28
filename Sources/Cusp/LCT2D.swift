//
//  LCT2D.swift
//  Cusp
//
//  Created by June Russell on 28/02/2025.
//

import Foundation
import Numerics
import Plinth

extension Matrix where Scalar == Double {
    
    public func lct2D(matrix: ComplexMatrix<Scalar>, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).lct2D(matrix: matrix, setup: setup)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func lct2D(matrix: ComplexMatrix<Scalar>, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        precondition(matrix.shape == .square(length: 2), "Transform matrix must be 2x2")
        
        var result = self
        
        // Apply LCT to each row (X dimension)
        for i in 0..<shape.rows {
            result[row: i] = result[row: i].lct1D(matrix: matrix, setup: setup)
        }
        
        // Apply LCT to each column (Y dimension)
        for j in 0..<shape.columns {
            result[column: j] = result[column: j].lct1D(matrix: matrix, setup: setup)
        }
        
        return result
    }
    
}