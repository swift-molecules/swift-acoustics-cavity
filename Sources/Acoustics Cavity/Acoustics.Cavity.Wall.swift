public import Mass

extension Acoustics.Cavity {
    /// A part of a cavity's boundary that moves, such as a soundboard or a bowl, given by its own modes in vacuo: the
    /// faces it covers and, for each mode, how each of their corners moves.
    ///
    /// A wall is whatever structure its modes come from, a plate, a shell or a membrane; the air needs only where
    /// it moves and how much mass and stiffness each mode has.
    public struct Wall {
        /// The faces the wall covers, each as three or four of the cavity's vertices, in the order of
        /// swift-galerkin's `Galerkin.Reference.Triangle` or `Quadrilateral`, turning counterclockwise seen from
        /// outside the air, so that their normal points out of it.
        public let faces: [[Int]]
        /// The wall's modes in vacuo.
        public let modes: [Mode]

        @inlinable
        public init(faces: [[Int]], modes: [Mode]) {
            self.faces = faces
            self.modes = modes
        }
    }
}

extension Acoustics.Cavity.Wall {
    /// One mode of a wall in vacuo, as the air sees it.
    public struct Mode {
        /// The natural frequency in vacuo.
        public let frequency: Frequency.Value<Double>
        /// The modal mass at the scale of `displacement`: the kinetic energy is one half this mass times the square of
        /// the rate of the mode's amplitude.
        public let mass: Mass.Value<Double>
        /// The loss factor η of the wall's own material, so that its stiffness is ω² m (1 + iη). Radiation is not
        /// in it: the air adds that.
        public let loss: Double
        /// The motion of each of the cavity's vertices, in its order, per unit of the mode's amplitude. Only the
        /// vertices of the wall's faces are read, and only the part along each face's normal moves the air, so a
        /// vertex where two faces meet at an angle, such as on a lute's rim, moves each as it should.
        public let displacement: [Vector<3, Double>]

        @inlinable
        public init(
            frequency: Frequency.Value<Double>,
            mass: Mass.Value<Double>,
            loss: Double,
            displacement: [Vector<3, Double>]
        ) {
            self.frequency = frequency
            self.mass = mass
            self.loss = loss
            self.displacement = displacement
        }
    }
}

extension Acoustics.Cavity.Wall: Hashable {}
extension Acoustics.Cavity.Wall: Sendable {}
extension Acoustics.Cavity.Wall.Mode: Hashable {}
extension Acoustics.Cavity.Wall.Mode: Sendable {}
