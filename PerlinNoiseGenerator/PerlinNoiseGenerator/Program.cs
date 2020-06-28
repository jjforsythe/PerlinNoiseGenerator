//Created by Max Whitby

//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "PerlinNoise"), to deal
//in the PerlinNoise without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:

//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.


using System;
using Gtk;
using Glade;

namespace PerlinNoiseGenerator
{
    class MainClass
    {
        public static void Main(string[] args) 
        {
            //Application.Init();
            //MainWindow win = new MainWindow();
            //win.Show();
            //Application.Run();

            Application.Init();

            //Create the Window
            Window myWin = new Window("Perlin Noise Generator v0.1");
            myWin.Resize(200, 200);
            myWin.Add(new DrawingArea());

            //Create a label and put some text in it.
            Label myLabel = new Label();
            myLabel.Text = "Hello World!!!!";

            //Add the label to the form
            myWin.Add(myLabel);

            //Show Everything
            myWin.ShowAll();

            Application.Run();
        }




    }
}
