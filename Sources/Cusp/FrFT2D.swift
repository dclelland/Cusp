//
//  FrFT2D.swift
//  Cusp
//
//  Created by June Russell on 28/02/2025.
//

import Foundation
import Numerics
import Plinth

extension Matrix where Scalar == Double {
    
    public func frft2D(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).frft2D(a: a, setup: setup)
    }
    
    public func frft2D(a_x: Scalar, a_y: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).frft2D(a_x: a_x, a_y: a_y, setup: setup)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    // Isotropic 2D FrFT (same angle in both dimensions)
    public func frft2D(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        return frft2D(a_x: a, a_y: a, setup: setup)
    }
    
    // Anisotropic 2D FrFT (different angles for each dimension)
    public func frft2D(a_x: Scalar, a_y: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        // First apply FrFT to each row
        let rowsTransformed = applyFrFTToRows(a: a_x, setup: setup)
        
        // Then apply FrFT to each column of the result
        return rowsTransformed.applyFrFTToColumns(a: a_y, setup: setup)
    }
    
    // Helper method to apply 1D FrFT to each row
    private func applyFrFTToRows(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        var result = ComplexMatrix<Scalar>.zeros(shape: shape)
        
        for i in 0..<shape.rows {
            let row = self[row: i]
            let transformedRow = row.frft1D(a: a, setup: setup)
            result[row: i] = transformedRow
        }
        
        return result
    }
    
    // Helper method to apply 1D FrFT to each column
    private func applyFrFTToColumns(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        var result = ComplexMatrix<Scalar>.zeros(shape: shape)
        
        for j in 0..<shape.columns {
            let column = self[column: j]
            let transformedColumn = column.frft1D(a: a, setup: setup)
            result[column: j] = transformedColumn
        }
        
        return result
    }
}
