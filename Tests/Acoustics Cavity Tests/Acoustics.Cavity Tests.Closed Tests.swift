import Acoustics_Cavity
import Testing

extension `Acoustics.Cavity Tests` {
    @Suite struct Closed {}
}

extension `Acoustics.Cavity Tests`.Closed {
    typealias Base = `Acoustics.Cavity Tests`

    /// f = (c / 2) √((l / x)² + (m / y)² + (n / z)²), every (l, m, n), sorted.
    static func box(_ x: Double, _ y: Double, _ z: Double, count: Int) -> [Double] {
        var result: [Double] = []
        for l in 0...6 {
            for m in 0...6 {
                for n in 0...6 {
                    let a = Double(l) / x
                    let b = Double(m) / y
                    let d = Double(n) / z
                    result.append(Base.c / 2 * (a * a + b * b + d * d).squareRoot())
                }
            }
        }
        return Array(result.sorted().prefix(count))
    }

    @Test
    func `A rigid box rings at its closed-form frequencies`() throws {
        let cavity = Base.Solid.box(0.5, 0.4, 0.3, cells: (20, 16, 12))
        let modes = try Acoustics.Cavity.Modes(cavity, count: 10, air: Base.air)
        let exact = Self.box(0.5, 0.4, 0.3, count: 10)
        #expect(modes[0].frequency.underlying < 1e-3)
        for (mode, f) in zip(modes.dropFirst(), exact.dropFirst()) {
            #expect(Base.error(mode.frequency.underlying, f) < 0.005, "\(mode.frequency.underlying) Hz, exact \(f)")
        }
    }

    @Test
    func `A rigid box of tetrahedra rings at its closed-form frequencies`() throws {
        let cavity = Base.Solid.box(0.5, 0.4, 0.3, cells: (15, 12, 9), tetrahedra: true)
        let modes = try Acoustics.Cavity.Modes(cavity, count: 8, air: Base.air)
        let exact = Self.box(0.5, 0.4, 0.3, count: 8)
        for (mode, f) in zip(modes.dropFirst(), exact.dropFirst()) {
            #expect(Base.error(mode.frequency.underlying, f) < 0.02, "\(mode.frequency.underlying) Hz, exact \(f)")
        }
    }

    @Test
    func `The box's frequencies converge at second order`() throws {
        // The (1, 1, 1) mode, refined twice.
        let exact = Self.box(0.5, 0.4, 0.3, count: 12)
        var errors: [Double] = []
        for n in [4, 8] {
            let cavity = Base.Solid.box(0.5, 0.4, 0.3, cells: (n * 5 / 4, n, n * 3 / 4))
            let modes = try Acoustics.Cavity.Modes(cavity, count: 12, air: Base.air)
            errors.append(Base.error(modes[11].frequency.underlying, exact[11]))
        }
        let order = (errors[0] / errors[1]).squareRoot()
        #expect(order > 1.8 && order < 2.2, "ratio \(errors[0] / errors[1])")
    }

    @Test
    func `A rigid cylinder rings at its closed-form frequencies`() throws {
        // f = (c / 2π) √((j′ / a)² + (p π / h)²), with j′ the zeros of Jₘ′: 1.841184 (m = 1), 3.054237 (m = 2),
        // 3.831706 (m = 0), each m > 0 twice.
        let a = 0.15
        let h = 0.3
        let cavity = Base.Solid.cylinder(radius: a, height: h, core: 5, rings: 5, layers: 16)
        let modes = try Acoustics.Cavity.Modes(cavity, count: 9, air: Base.air)
        func f(_ j: Double, _ p: Int) -> Double {
            let axial = Double(p) * Double.pi / h
            return Base.c / (2 * Double.pi) * ((j / a) * (j / a) + axial * axial).squareRoot()
        }
        let exact = [
            0, f(0, 1), f(1.841184, 0), f(1.841184, 0), f(1.841184, 1), f(1.841184, 1), f(0, 2), f(3.054237, 0),
            f(3.054237, 0),
        ]
        .sorted()
        for (mode, f) in zip(modes.dropFirst(), exact.dropFirst()) {
            #expect(Base.error(mode.frequency.underlying, f) < 0.01, "\(mode.frequency.underlying) Hz, exact \(f)")
        }
    }

    @Test
    func `A rigid sphere rings at its closed-form frequencies`() throws {
        // f = c x / 2π a, with x the zeros of jₗ′: 2.081576 (l = 1, three modes), 3.342094 (l = 2, five),
        // 4.493409 (l = 0, one), 4.514100 (l = 3, seven).
        let a = 0.2
        let cavity = Base.Solid.ball(radius: a, core: 4, shells: 4)
        let modes = try Acoustics.Cavity.Modes(cavity, count: 17, air: Base.air)
        let roots =
            [2.081576, 2.081576, 2.081576] + Array(repeating: 3.342094, count: 5) + [4.493409]
            + Array(repeating: 4.514100, count: 7)
        #expect(modes[0].frequency.underlying < 1e-3)
        for (mode, x) in zip(modes.dropFirst(), roots) {
            let f = Base.c * x / (2 * Double.pi * a)
            #expect(Base.error(mode.frequency.underlying, f) < 0.01, "\(mode.frequency.underlying) Hz, exact \(f)")
        }
    }
}
