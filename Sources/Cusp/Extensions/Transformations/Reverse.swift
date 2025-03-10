//
//  Reverse.swift
//  Cusp
//
//  Created by Daniel Clelland on 07/03/2025.
//

import Foundation
import Plinth

extension ComplexMatrix where Scalar == Float {
    
    internal func reversedBody() -> ComplexMatrix {
        return ComplexMatrix(real: real.reversedBody(), imaginary: imaginary.reversedBody())
    }
    
}

extension Matrix where Scalar == Float {
    
    internal func reversedBody() -> Matrix {
        let head = elements.first!
        let body = elements.dropFirst()
        return Matrix(shape: shape, elements: [head] + body.reversed())
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    internal func reversedBody() -> ComplexMatrix {
        return ComplexMatrix(real: real.reversedBody(), imaginary: imaginary.reversedBody())
    }
    
}

extension Matrix where Scalar == Double {
    
    internal func reversedBody() -> Matrix {
        let head = elements.first!
        let body = elements.dropFirst()
        return Matrix(shape: shape, elements: [head] + body.reversed())
    }
    
}

extension Matrix where Scalar == Double {
    
    public func reversedRows() -> Matrix {
        var output = self
        for row in output.shape.rowIndices {
            output[row: row] = output[row: row].reversed()
        }
        return output
    }
    
    public func reversedColumns() -> Matrix {
        var output = self
        for column in output.shape.columnIndices {
            output[column: column] = output[column: column].reversed()
        }
        return output
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func reversedRows() -> ComplexMatrix {
        return ComplexMatrix(real: real.reversedRows(), imaginary: imaginary.reversedRows())
    }
    
    public func reversedColumns() -> ComplexMatrix {
        return ComplexMatrix(real: real.reversedColumns(), imaginary: imaginary.reversedColumns())
    }
    
}
