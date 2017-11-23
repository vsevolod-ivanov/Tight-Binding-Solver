

//C Program to Illustrate Reading of Data from a File

//This C Program illustrates reading of data from a file. The program opens a file which is present. Once the file opens successfully, it uses libc fgetc() library call to read the content.
//Here is source code of the C program to illustrate reading of data from a file. The C program is successfully compiled and run on a Linux system. The program output is also shown below.

/*
 * C program to illustrate how a file stored on the disk is read
 */
#include <stdio.h>
#include <stdlib.h>
 
void main()
{
    FILE *fptr;
    char filename[15];
    char ch;
 
    printf("Enter the filename to be opened \n");
    scanf("%s", filename);
    /*  open the file for reading */
    fptr = fopen(filename, "r");
    if (fptr == NULL)
    {
        printf("Cannot open file \n");
        exit(0);
    }
    ch = fgetc(fptr);
    while (ch != EOF)
    {
        printf ("%c", ch);
        ch = fgetc(fptr);
    }
    fclose(fptr);
}
//$ cc pgm96.c
//$ a.out
//Enter the filename to be opened
//pgm95.c
/*
 * C program to create a file called emp.rec and store information
 * about a person, in terms of his name, age and salary.
 */
 
#include <stdio.h>
 
void main()
{
    FILE *fptr;
    char name[20];
    int age;
    float salary;
 
    fptr = fopen ("emp.rec", "w"); /* open for writing*/
 
    if (fptr == NULL)
    {
        printf("File does not exists \n");
        return;
    }
    printf("Enter the name \n");
    scanf("%s", name);
    fprintf(fptr, "Name    = %s\n", name);
    printf("Enter the age \n");
    scanf("%d", &age);
    fprintf(fptr, "Age     = %d\n", age);
    printf("Enter the salary \n");
    scanf("%f", &salary);
    fprintf(fptr, "Salary  = %.2f\n", salary);
    fclose(fptr);
}