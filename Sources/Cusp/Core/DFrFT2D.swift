//
//  DFrFT2D.swift
//  Cusp
//
//  Created by Daniel Clelland on 09/03/2025.
//

import Foundation
import Plinth

extension Matrix where Scalar == Double {
    
    public func dfrft2D(order: Scalar, approximationOrder: Int = 2) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).dfrft2D(order: order, approximationOrder: approximationOrder)
    }

    public func dfrft2D(matrix: ComplexMatrix<Scalar>) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).dfrft2D(matrix: matrix)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func dfrft2D(order: Scalar, approximationOrder: Int = 2) -> ComplexMatrix<Scalar> {
        let matrix = ComplexMatrix.dfrftMatrix(length: shape.count, order: order, approximationOrder: approximationOrder)
        return dfrft2D(matrix: matrix)
    }
    
    public func dfrft2D(matrix: ComplexMatrix<Scalar>) -> ComplexMatrix<Scalar> {
        var result = self
        
        for i in 0..<shape.rows {
            result[row: i] = result[row: i].dfrft1D(matrix: matrix)
        }
        
        for j in 0..<shape.columns {
            result[column: j] = result[column: j].asRow().dfrft1D(matrix: matrix).asColumn()
        }
        
        return result
    }
    
}
