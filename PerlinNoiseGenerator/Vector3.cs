﻿using System;
using System.Collections.Generic;

namespace max.linearalgebra
{
    public class Vector3<T>
    {
        protected T[] coord;

        public Vector3()
        {
            coord = new T[3];
        }

        public Vector3(T x1, T x2, T x3)
        {
            coord = new T[3] { x1, x2, x3 };
        }

        public Vector3(Vector3<T> v)
        {
            coord = new T[3] { v[0], v[1], v[2] };
        }

        public Vector3(T[] v)
        {
            coord = new T[3] { v[0], v[1], v[2] };
        }

        public T this[int i]
        {
            get
            {
                return coord[i];
            }
            set
            {
                coord[i] = value;
            }
        }


        public T X
        {
            get
            {
                return coord[0];
            }
            set
            {
                coord[0] = value;
            }
        }

        public T Y
        {
            get
            {
                return coord[1];
            }
            set
            {
                coord[1] = value;
            }
        }

        public T Z
        {
            get
            {
                return coord[2];
            }
            set
            {
                coord[2] = value;
            }
        }

    }

    public class PhysicalVector3 : vector3<double>
    {
        public PhysicalVector3(double x1, double x2, double x3) : base(x1, x2, x3) { }
        public PhysicalVector3(PhysicalVector3 v) : base(v) { }
        public PhysicalVector3(double[] v) : base(v) { }
        public PhysicalVector3() : base() { }


        // l2-norm for vectors = sqrt(x1^2+x2^2+x3^2)
        public double L2norm() => Math.Sqrt(coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2]);


        // STOP GAP PLEASE REFACTOR LOWER CASE VERSION OUT AT SOME POINT
        public double l2norm() => Math.Sqrt(coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2]);

        // square of the l2-norm for vectors = (x1^2+x2^2+x3^2) - sqrt is expensive
        public double L2norm2() => coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2];



        // Needs serious refactor.. why does Normalize return and Invert doesnt???
        // Inserted after to support mesh simplification
        public PhysicalVector3 Normalize()
        {
            // Implements magnitude function
            double mag = L2norm();

            // Normalize with respect to the magnitude
            PhysicalVector3 result = new PhysicalVector3();

            if (mag > double.Epsilon)
            {
                result = new PhysicalVector3(coord[0] / mag, coord[1] / mag, coord[2] / mag);
            }
            else
            {
                result = new PhysicalVector3(0, 0, 0);
            }
            return result;
        }

        public void Invert()
        {
            this.X = -X; this.Y = -Y; this.Z = -Z;
        }

        public static double Dot(PhysicalVector3 v1, PhysicalVector3 v2) => (v1.X * v2.X) + (v1.Y * v2.Y) + (v1.Z * v2.Z);


        public static PhysicalVector3 Tangent(PhysicalVector3 v1, PhysicalVector3 v2) => new PhysicalVector3(v2.X - v1.X, v2.Y - v1.Y, v2.Z - v1.Z);


        public static PhysicalVector3 Cross(PhysicalVector3 v1, PhysicalVector3 v2) => new PhysicalVector3(v1.Y * v2.Z - v1.Z * v2.Y, v1.Z * v2.X - v1.X * v2.Z, v1.X * v2.Y - v1.Y * v2.X);

        public static bool ApproxEqual(PhysicalVector3 v1, PhysicalVector3 v2, double tolerance) => (Math.Abs(v1.X - v2.X) < tolerance && Math.Abs(v1.Y - v2.Y) < tolerance
                                                                                        && Math.Abs(v1.Z - v2.Z) < tolerance) ? true : false;


        // l2-norm-based squared distance between vectors
        public static double Dist(PhysicalVector3 v1, PhysicalVector3 v2) => Math.Sqrt(Dist2(v1, v2));


        // square of the l2-norm-based squared distance between vectors - sqrt is expensive
        public static double Dist2(PhysicalVector3 v1, PhysicalVector3 v2) => (v1[0] - v2[0]) * (v1[0] - v2[0])
                                                         + (v1[1] - v2[1]) * (v1[1] - v2[1])
                                                         + (v1[2] - v2[2]) * (v1[2] - v2[2]);


        // scalar product of two 3-vectors
        public static double ScalarProduct(PhysicalVector3 v1, PhysicalVector3 v2) => v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];


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
            double nn = nml.L2norm();
            double n0 = nml[0] / nn, n1 = nml[1] / nn, n2 = nml[2] / nn;
            double scaling = 1.0 / (1.0 + n2);

            xi.Set(1 - n0 * n0 * scaling, -n0 * n1 * scaling, -n0);
            eta.Set(-n0 * n1 * scaling, 1 - n1 * n1 * scaling, -n1);
        }

        public static void CrossAxesX(PhysicalVector3 nml, ref PhysicalVector3 xi, ref PhysicalVector3 eta)
        {
            double nn = nml.L2norm();
            double n0 = nml[0] / nn, n1 = nml[1] / nn, n2 = nml[2] / nn;
            double scaling = 1.0 / (1.0 + n0);

            xi.Set(-n1, 1 - n1 * n1 * scaling, -n1 * n2 * scaling);
            eta.Set(-n2, -n1 * n2 * scaling, 1 - n2 * n2 * scaling);
        }


        public static double Angle(PhysicalVector3 v1, PhysicalVector3 v2)
        {
            PhysicalVector3 fromNormalized = v1.Normalize();
            PhysicalVector3 toNormalized = v2.Normalize();
            return Math.Acos(MathUtils.Clamp(PhysicalVector3.ScalarProduct(fromNormalized, toNormalized), -1.0, 1.0)) * MathUtils.Rad2Degd;
        }

        public static double signedAngle(PhysicalVector3 v1, PhysicalVector3 v2)
        {


            return Math.Atan2(PhysicalVector3.Cross(v1, v2).L2norm2(), PhysicalVector3.Dot(v1, v2));

        }


        // Rodrigues Formula to rotate vector V about Axis of angle theta
        public static void RotateVector(PhysicalVector3 V, PhysicalVector3 Axis, double theta, ref PhysicalVector3 Vrot)
        {
            double cosTheta = Math.Cos(theta);
            double sinTheta = Math.Sin(theta);
            PhysicalVector3 U1 = new PhysicalVector3(V[0] * cosTheta, V[1] * cosTheta, V[2] * cosTheta);
            CrossProduct(Axis, V, out PhysicalVector3 U2);
            U2.Set(U2[0] * sinTheta, U2[1] * sinTheta, U2[2] * sinTheta);
            double dot = ScalarProduct(Axis, V);
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

        public static PhysicalVector3 operator /(PhysicalVector3 vector1, double d1)
        {
            return new PhysicalVector3(vector1.X / d1, vector1.Y / d1, vector1.Z / d1);
        }

        public static PhysicalVector3 operator *(PhysicalVector3 vector1, double d1)
        {
            return new PhysicalVector3(vector1.X * d1, vector1.Y * d1, vector1.Z * d1);
        }

        public override string ToString() { return "{ X: " + X.ToString() + " Y: " + Y.ToString() + " Z: " + Z.ToString() + "} "; }




    }


    //
    // dvec3-int tuple
    //
    public struct dvec3int
    {
        public double x;
        public double y;
        public double z;
        public int nodeind;

        public dvec3int(PhysicalVector3 point, int ni)
        {
            x = point[0];
            y = point[1];
            z = point[2];
            nodeind = ni;
        }

        public void Set(PhysicalVector3 point, int ni)
        {
            x = point[0];
            y = point[1];
            z = point[2];
            nodeind = ni;
        }
    }

    //
    // Dictionary with key = index of an edge in the 3D image
    //                 value = coordinate of the point where this edge is intersected by phi
    //
    public class IntToDvec3IntDictionary : Dictionary<int, dvec3int>
    {
        public readonly int[] dims;

        public IntToDvec3IntDictionary(int nx, int ny, int nz, int capacity)
            : base(capacity)
        {
            dims = new int[] { nx, ny, nz };
        }
    }


    //
    // Class used for comparison of double 3-vectors
    // Distance between grid vertices is O(1), 
    // so rounding up to 2 digits should be sufficient to identify intersection points uniquely
    //
    public class NodeEqualityComparer : IEqualityComparer<PhysicalVector3>
    {

        // Nodes are combined if distance between them is less than 1e-2
        public bool Equals(PhysicalVector3 v1, PhysicalVector3 v2)
        {
            return (PhysicalVector3.Dist2(v1, v2) < MathUtils.DOUBLE_TOLERANCE);
        }

        // Hash by nearest integer to the norm 
        // In image space node coordinates are between 0 and O(100), so precision of O(1) should be sufficient
        public int GetHashCode(PhysicalVector3 v)
        {
            return (int)v.L2norm2();
        }
    }



}
