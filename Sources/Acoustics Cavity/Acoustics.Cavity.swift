extension Acoustics {
    /// A body of air: a volume cut into hexahedra or tetrahedra, held by rigid walls except where its pressure is
    /// released.
    ///
    /// A closed body, such as a box, has rigid walls everywhere. A body with an opening, such as a guitar with its
    /// sound hole or a lute with its rose, is given with some of the air outside the opening, and that outer air's
    /// far boundary is `open`.
    public struct Cavity {
        /// The vertices, in metres.
        public let vertices: [Point<3, Double>]
        /// Each cell's corners as indices into `vertices`: eight for a hexahedron, in the order of swift-galerkin's
        /// `Galerkin.Reference.Hexahedron`, or four for a tetrahedron, in the order of its `Tetrahedron`. A cell
        /// whose corners turn the other way is accepted and turned.
        public let cells: [[Int]]
        /// The vertices where the pressure is released to zero; everywhere else the boundary is a rigid wall.
        public let open: Set<Int>

        @inlinable
        public init(vertices: [Point<3, Double>], cells: [[Int]], open: Set<Int> = []) {
            self.vertices = vertices
            self.cells = cells
            self.open = open
        }
    }
}

extension Acoustics.Cavity: Hashable {}
extension Acoustics.Cavity: Sendable {}
