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
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func frft2D(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        let rowsTransformed = applyFrFTToRows(a: a, setup: setup)
        return rowsTransformed.applyFrFTToColumns(a: a, setup: setup)
    }
    
    private func applyFrFTToRows(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        var result = ComplexMatrix<Scalar>.zeros(shape: shape)
        
        for i in 0..<shape.rows {
            let row = self[row: i]
            let transformedRow = row.frft1D(a: a, setup: setup)
            result[row: i] = transformedRow
        }
        
        return result
    }
    
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
