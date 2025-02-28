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
    
    public static func fractionalFourierTransform(a: Scalar) -> ComplexMatrix<Scalar> {
        let alpha = a * .pi / 2.0
        return [
            [.init(cos(alpha)), .init(sin(alpha))],
            [.init(-sin(alpha)), .init(cos(alpha))]
        ]
    }
    
}
