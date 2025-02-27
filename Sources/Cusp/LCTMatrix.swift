//
//  LCTMatrix.swift
//  Cusp
//
//  Created by June Russell on 28/02/2025.
//

import Foundation
import Numerics
import Plinth

public enum LCTMatrix<Scalar> { }

extension LCTMatrix where Scalar == Double {
    
    public static let identityTransform: ComplexMatrix<Scalar> = [
        [.one, .zero],
        [.zero, .one]
    ]
    
    public static let fourierTransform: ComplexMatrix<Scalar> = [
        [.zero, .one],
        [-.one, .zero]
    ]
    
    public static let inverseFourierTransform: ComplexMatrix<Scalar> = [
        [.zero, -.one],
        [.one, .zero]
    ]
    
    public static func fractionalFourierTransform(angle: Scalar) -> ComplexMatrix<Scalar> {
        return [
            [.init(cos(angle)), .init(sin(angle))],
            [.init(-sin(angle)), .init(cos(angle))]
        ]
    }
    
}
