import Acoustics_Cavity
import Testing

extension `Acoustics.Cavity Tests` {
    @Suite struct Laws {}
}

extension `Acoustics.Cavity Tests`.Laws {
    typealias Base = `Acoustics.Cavity Tests`

    static let box = Base.Solid.box(0.5, 0.4, 0.3, cells: (10, 8, 6))

    @Test
    func `A closed cavity's first mode is a uniform pressure whose modal volume is the cavity's`() throws {
        let modes = try Acoustics.Cavity.Modes(Self.box, count: 2, air: Base.air)
        #expect(modes[0].frequency.underlying < 1e-3)
        #expect(modes[0].shape.allSatisfy { Swift.abs($0 - 1) < 1e-9 })
        #expect(Swift.abs(modes[0].volume.underlying - 0.5 * 0.4 * 0.3) < 1e-12)
        // The (1, 0, 0) mode is cos(π x / 0.5): its modal volume is half the box's, less the mesh's error of
        // about (k h)² / 6 with ten cells along x.
        #expect(Base.error(modes[1].volume.underlying, 0.03) < 0.02)
    }

    @Test
    func `Modes ascend and each shape's largest pressure is one`() throws {
        let modes = try Acoustics.Cavity.Modes(Self.box, count: 8, air: Base.air)
        #expect(modes.count == 8)
        for (lower, upper) in zip(modes, modes.dropFirst()) {
            #expect(lower.frequency.underlying <= upper.frequency.underlying)
        }
        for mode in modes {
            #expect(mode.shape.count == Self.box.vertices.count)
            #expect(mode.shape.max() == 1)
            #expect(mode.shape.allSatisfy { $0 >= -1 && $0 <= 1 })
        }
    }

    @Test
    func `Open vertices have no pressure`() throws {
        // Release the face x = 0: the box becomes a closed–open duct along x, f = c (2l + 1) / 4 x at the lowest.
        let closed = Self.box
        let open = Set(closed.vertices.indices.filter { closed.vertices[$0].x == 0 })
        let cavity = Acoustics.Cavity(vertices: closed.vertices, cells: closed.cells, open: open)
        let modes = try Acoustics.Cavity.Modes(cavity, count: 3, air: Base.air)
        for mode in modes { #expect(open.allSatisfy { mode.shape[$0] == 0 }) }
        #expect(Base.error(modes[0].frequency.underlying, Base.c / (4 * 0.5)) < 0.005)
    }

    @Test
    func `A cavity scaled by s rings at 1 / s the frequency with s³ the modal volume`() throws {
        let s = 2.5
        let scaled = Acoustics.Cavity(
            vertices: Self.box.vertices.map { Point(x: s * $0.x, y: s * $0.y, z: s * $0.z) },
            cells: Self.box.cells
        )
        let modes = try Acoustics.Cavity.Modes(Self.box, count: 6, air: Base.air)
        let larger = try Acoustics.Cavity.Modes(scaled, count: 6, air: Base.air)
        for (mode, other) in zip(modes.dropFirst(), larger.dropFirst()) {
            #expect(Base.error(other.frequency.underlying * s, mode.frequency.underlying) < 1e-9)
            #expect(Base.error(other.volume.underlying, s * s * s * mode.volume.underlying) < 1e-9)
        }
    }

    @Test
    func `Frequencies are proportional to the speed of sound`() throws {
        var warm = Base.air
        warm.speed = Velocity.amount(1.1 * Base.c)
        let modes = try Acoustics.Cavity.Modes(Self.box, count: 4, air: Base.air)
        let faster = try Acoustics.Cavity.Modes(Self.box, count: 4, air: warm)
        for (mode, other) in zip(modes.dropFirst(), faster.dropFirst()) {
            #expect(Base.error(other.frequency.underlying, 1.1 * mode.frequency.underlying) < 1e-12)
            #expect(other.shape == mode.shape && other.volume == mode.volume)
        }
    }

    @Test
    func `A cell listed the other way round is turned`() throws {
        let hexahedra = Self.box.cells.map { cell in [4, 5, 6, 7, 0, 1, 2, 3].map { cell[$0] } }
        let turned = Acoustics.Cavity(vertices: Self.box.vertices, cells: hexahedra)
        let modes = try Acoustics.Cavity.Modes(Self.box, count: 5, air: Base.air)
        let other = try Acoustics.Cavity.Modes(turned, count: 5, air: Base.air)
        for (mode, again) in zip(modes.dropFirst(), other.dropFirst()) {
            #expect(Base.error(again.frequency.underlying, mode.frequency.underlying) < 1e-12)
        }
        let tetrahedra = Base.Solid.box(0.5, 0.4, 0.3, cells: (5, 4, 3), tetrahedra: true)
        let swapped = Acoustics.Cavity(
            vertices: tetrahedra.vertices,
            cells: tetrahedra.cells.map { [$0[0], $0[2], $0[1], $0[3]] }
        )
        let lower = try Acoustics.Cavity.Modes(tetrahedra, count: 4, air: Base.air)
        let again = try Acoustics.Cavity.Modes(swapped, count: 4, air: Base.air)
        for (mode, other) in zip(lower.dropFirst(), again.dropFirst()) {
            #expect(Base.error(other.frequency.underlying, mode.frequency.underlying) < 1e-12)
        }
    }

    @Test
    func `Bad cells and counts fail`() {
        let vertices = Self.box.vertices
        let good = Self.box.cells[0]
        #expect(throws: Acoustics.Cavity.Error.cell(1)) {
            try Acoustics.Cavity.Modes(Acoustics.Cavity(vertices: vertices, cells: [good, [0, 1, 2]]), count: 1)
        }
        #expect(throws: Acoustics.Cavity.Error.cell(0)) {
            try Acoustics.Cavity.Modes(
                Acoustics.Cavity(vertices: vertices, cells: [[0, 1, 2, vertices.count]]), count: 1)
        }
        // A tetrahedron with no volume, and a hexahedron folded through itself.
        #expect(throws: Acoustics.Cavity.Error.cell(0)) {
            try Acoustics.Cavity.Modes(Acoustics.Cavity(vertices: vertices, cells: [[0, 1, 2, 3]]), count: 1)
        }
        let folded = [good[0], good[1], good[3], good[2], good[4], good[5], good[6], good[7]]
        #expect(throws: Acoustics.Cavity.Error.cell(0)) {
            try Acoustics.Cavity.Modes(Acoustics.Cavity(vertices: vertices, cells: [folded]), count: 1)
        }
        #expect(throws: Acoustics.Cavity.Error.count) {
            try Acoustics.Cavity.Modes(Self.box, count: Self.box.vertices.count + 1)
        }
    }
}
