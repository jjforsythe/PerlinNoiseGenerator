using System;
namespace max.linearalgebra
{

    public class PhysicalVector3 : Vector3<float>
    {
        public PhysicalVector3(float x1, float x2, float x3) : base(x1, x2, x3) { }
        public PhysicalVector3(PhysicalVector3 v) : base(v) { }
        public PhysicalVector3(float[] v) : base(v) { }
        public PhysicalVector3() : base() { }


        // l2-norm for vectors = sqrt(x1^2+x2^2+x3^2)
        public float L2norm() => Math.Sqrt(coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2]);

        // square of the l2-norm for vectors = (x1^2+x2^2+x3^2) - sqrt is expensive
        public float L2norm2() => coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2];


        public void Normalize()
        {
            // Implements magnitude function
            float mag = L2norm();

            // Normalize with respect to the magnitude
            PhysicalVector3 result = new PhysicalVector3();

            if (mag > float.Epsilon)
            {
                result = new PhysicalVector3(coord[0] / mag, coord[1] / mag, coord[2] / mag);
            }
            else
            {
                result = new PhysicalVector3(0, 0, 0);
            }
        }

        public void Invert()
        {
            this.X = -X; this.Y = -Y; this.Z = -Z;
        }

        public static float Dot(PhysicalVector3 v1, PhysicalVector3 v2) => (v1.X * v2.X) + (v1.Y * v2.Y) + (v1.Z * v2.Z);


        public static PhysicalVector3 Tangent(PhysicalVector3 v1, PhysicalVector3 v2) => new PhysicalVector3(v2.X - v1.X, v2.Y - v1.Y, v2.Z - v1.Z);


        public static PhysicalVector3 Cross(PhysicalVector3 v1, PhysicalVector3 v2) => new PhysicalVector3(v1.Y * v2.Z - v1.Z * v2.Y, v1.Z * v2.X - v1.X * v2.Z, v1.X * v2.Y - v1.Y * v2.X);

        public static bool ApproxEqual(PhysicalVector3 v1, PhysicalVector3 v2, float tolerance) => (Math.Abs(v1.X - v2.X) < tolerance && Math.Abs(v1.Y - v2.Y) < tolerance
                                                                                        && Math.Abs(v1.Z - v2.Z) < tolerance) ? true : false;


        // l2-norm-based squared distance between vectors
        public static float Dist(PhysicalVector3 v1, PhysicalVector3 v2) => Math.Sqrt(Dist2(v1, v2));


        // square of the l2-norm-based squared distance between vectors - sqrt is expensive
        public static float Dist2(PhysicalVector3 v1, PhysicalVector3 v2) => (v1[0] - v2[0]) * (v1[0] - v2[0])
                                                         + (v1[1] - v2[1]) * (v1[1] - v2[1])
                                                         + (v1[2] - v2[2]) * (v1[2] - v2[2]);


        // scalar product of two 3-vectors
        public static float ScalarProduct(PhysicalVector3 v1, PhysicalVector3 v2) => v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];


        // Inserted to support mesh simplification 
        //cross product
        //So bad... Please 
        public static void CrossProduct(PhysicalVector3 v1, PhysicalVector3 v2, out PhysicalVector3 result)
        {
            //lhs = v1, rhs = v2
            //lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x
            result = new PhysicalVector3(v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]);
        }

        // Consistently find orthonormal unit vectors <xi,eta> 
        // in the plane defined by normal nml
        // by rotating <i,j,k> into <xi,eta,nml>
        public static void CrossAxesZ(PhysicalVector3 nml, ref PhysicalVector3 xi, ref PhysicalVector3 eta)
        {
            float nn = nml.L2norm();
            float n0 = nml[0] / nn, n1 = nml[1] / nn, n2 = nml[2] / nn;
            float scaling = 1.0 / (1.0 + n2);

            xi.Set(1 - n0 * n0 * scaling, -n0 * n1 * scaling, -n0);
            eta.Set(-n0 * n1 * scaling, 1 - n1 * n1 * scaling, -n1);
        }

        public static void CrossAxesX(PhysicalVector3 nml, ref PhysicalVector3 xi, ref PhysicalVector3 eta)
        {
            float nn = nml.L2norm();
            float n0 = nml[0] / nn, n1 = nml[1] / nn, n2 = nml[2] / nn;
            float scaling = 1.0 / (1.0 + n0);

            xi.Set(-n1, 1 - n1 * n1 * scaling, -n1 * n2 * scaling);
            eta.Set(-n2, -n1 * n2 * scaling, 1 - n2 * n2 * scaling);
        }


        public static float Angle(PhysicalVector3 v1, PhysicalVector3 v2)
        {
            PhysicalVector3 fromNormalized = v1.Normalize();
            PhysicalVector3 toNormalized = v2.Normalize();
            return Math.Acos(MathUtils.Clamp(PhysicalVector3.ScalarProduct(fromNormalized, toNormalized), -1.0, 1.0)) * MathUtils.Rad2Degd;
        }

        public static float signedAngle(PhysicalVector3 v1, PhysicalVector3 v2)
        {


            return Math.Atan2(PhysicalVector3.Cross(v1, v2).L2norm2(), PhysicalVector3.Dot(v1, v2));

        }


        // Rodrigues Formula to rotate vector V about Axis of angle theta
        public static void RotateVector(PhysicalVector3 V, PhysicalVector3 Axis, float theta, ref PhysicalVector3 Vrot)
        {
            float cosTheta = Math.Cos(theta);
            float sinTheta = Math.Sin(theta);
            PhysicalVector3 U1 = new PhysicalVector3(V[0] * cosTheta, V[1] * cosTheta, V[2] * cosTheta);
            CrossProduct(Axis, V, out PhysicalVector3 U2);
            U2.Set(U2[0] * sinTheta, U2[1] * sinTheta, U2[2] * sinTheta);
            float dot = ScalarProduct(Axis, V);
            PhysicalVector3 U3 = new PhysicalVector3(Axis[0] * dot * (1 - cosTheta), Axis[1] * dot * (1 - cosTheta), Axis[2] * dot * (1 - cosTheta));
            Vrot = new PhysicalVector3(U1[0] + U2[0] + U3[0], U1[1] + U2[1] + U3[1], U1[2] + U2[2] + U3[2]);
        }



        public static PhysicalVector3 operator -(PhysicalVector3 vector1, PhysicalVector3 vector2)
        {
            return new PhysicalVector3(vector1.X - vector2.X, vector1.Y - vector2.Y, vector1.Z - vector2.Z);
        }

        public static PhysicalVector3 operator -(PhysicalVector3 vector1)
        {
            return new PhysicalVector3(-vector1.X, -vector1.Y, -vector1.Z);
        }

        public static PhysicalVector3 operator +(PhysicalVector3 vector1, PhysicalVector3 vector2)
        {
            return new PhysicalVector3(vector1.X + vector2.X, vector1.Y + vector2.Y, vector1.Z + vector2.Z);
        }

        public static PhysicalVector3 operator /(PhysicalVector3 vector1, float d1)
        {
            return new PhysicalVector3(vector1.X / d1, vector1.Y / d1, vector1.Z / d1);
        }

        public static PhysicalVector3 operator *(PhysicalVector3 vector1, float d1)
        {
            return new PhysicalVector3(vector1.X * d1, vector1.Y * d1, vector1.Z * d1);
        }

        public override string ToString() { return "{ X: " + X.ToString() + " Y: " + Y.ToString() + " Z: " + Z.ToString() + "} "; }


    }

}
