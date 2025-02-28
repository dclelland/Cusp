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
    
    public func lct2D(matrixX: ComplexMatrix<Scalar>, matrixY: ComplexMatrix<Scalar>, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).lct2D(matrixX: matrixX, matrixY: matrixY, setup: setup)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    // Isotropic 2D LCT (same matrix in both dimensions)
    public func lct2D(matrix: ComplexMatrix<Scalar>, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        return lct2D(matrixX: matrix, matrixY: matrix, setup: setup)
    }
    
    // Anisotropic 2D LCT (different matrices for each dimension)
    public func lct2D(matrixX: ComplexMatrix<Scalar>, matrixY: ComplexMatrix<Scalar>, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        precondition(matrixX.shape == .square(length: 2), "X transform matrix must be 2x2")
        precondition(matrixY.shape == .square(length: 2), "Y transform matrix must be 2x2")
        
        var result = self
        
        // Apply LCT to each row (X dimension)
        for i in 0..<shape.rows {
            result[row: i] = result[row: i].lct1D(matrix: matrixX, setup: setup)
        }
        
        // Apply LCT to each column (Y dimension)
        for j in 0..<shape.columns {
            result[column: j] = result[column: j].lct1D(matrix: matrixY, setup: setup)
        }
        
        return result
    }
    
}