using System;
namespace Application
{

    public static class ImageUtils
    {

        public static Vector3[,] IntensityToRGB(float[] image) => IntensityToRGB(image, image.GetLength(0), image.GetLength(1));

        public static Vector3[,] IntensityToRGB(float[] image, int width, int height)
        {
            Vector3[,] bitmap = new Vector3[width, height];

        }
    }
}
