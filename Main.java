import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;

public class Main {
    static double[][] Matrix;
    static int size;
    static double[][] U;
    static double[][] L;
    static double[][] P;
    static double[][] Ptilda;
    static ArrayList<rowsInPs> rowsInPs = new ArrayList<>();
    static ArrayList<double[][]> Ms = new ArrayList<>();
    static ArrayList<double[][]> Mtildas = new ArrayList<>();

    public static void main(String[] args) {
        long startTime = System.nanoTime();
        size = ThreadLocalRandom.current().nextInt(10, 300);
//        double[][] A =random(size,size);
        double[][] A =  subtract(multiply(random(size, size - 1), random(size - 1, size)), rI(-0.03, identity()));
        Matrix = new double[size][size];
        copy(A, Matrix);
        P = identity();
        Ptilda = identity();
        LU_Factorization(Matrix);
        long endTime = System.nanoTime();
        System.out.println("Size = " + size);
        NumberFormat formatter = new DecimalFormat("#0.000000000");
        System.out.println("Execution time is " + formatter.format((endTime - startTime) / 1000000000d) + " seconds");
        System.out.println("||PAP~ - LU||  = " + frobeniusNorm(subtract(multiply(multiply(P, A), Ptilda), multiply(L, U))));
    }

    static class rowsInPs {
        public int rowA;
        public int rowB;

        rowsInPs(int a, int b) {
            this.rowA = a;
            this.rowB = b;
        }
    }

    static double[][] random(int numberOfRowa, int numOfColumns) {
        double[][] a = new double[numberOfRowa][numOfColumns];
        int min = ThreadLocalRandom.current().nextInt(-1000, 0),
                max = ThreadLocalRandom.current().nextInt(0, 1000);
        System.out.println("min = " + min + "\t,  max = " + max);
        for (int i = 0; i < numberOfRowa; i++) {
            for (int j = 0; j < numOfColumns; j++) {
                a[i][j] = ThreadLocalRandom.current().nextDouble(min, max);

            }
        }
        return a;
    }

    static double[][] rI(double r, double[][] I) { // return r*I
        for (int i = 0; i < size; i++)
            I[i][i] = r;
        return I;
    }

    static void copy(double[][] a, double[][] b) {
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                b[i][j] = a[i][j];
    }

    static void LU_Factorization(double[][] mat) {
        for (int k = 0; k < size - 1; k++) {
            int rowOfPivot = k, columnOfPivot = k;
            for (int i = k; i < size; i++) { //find largest elements in matrix[k:n,k:n]
                for (int j = k; j < size; j++) {
                    if (Math.abs(mat[i][j]) > Math.abs(mat[rowOfPivot][columnOfPivot])) {
                        rowOfPivot = i;
                        columnOfPivot = j;
                    }
                }
            }
            rowsInPs.add(new rowsInPs(k, rowOfPivot));
            changeRows(P, rowOfPivot, k); // P = p_k * P
            changeColumns(Ptilda, columnOfPivot, k); // Ptilda = Ptilda * p_k
            changeRows(mat, rowOfPivot, k); //update matrix A
            changeColumns(mat, columnOfPivot, k);
            double largest = mat[k][k];
            double[][] M = identity();
            for (int t = k + 1; t < size; t++)
                M[t][k] = -mat[t][k] / largest;
            Ms.add(M);
            for (int t = k + 1; t < size; t++) //update matrix A
                for (int q = k + 1; q < size; q++)
                    mat[t][q] = mat[t][q] + M[t][k] * mat[k][q];
        }
        U = new double[size][size]; // create U
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++)
                if (j >= i)
                    U[i][j] = mat[i][j];
        }
        createL();
    }

    static void createL() {
        L = identity();
        int i = size - 2;
        double[][] pl = identity();
        changeRows(pl, rowsInPs.get(i).rowA, rowsInPs.get(i).rowB);
        double[][] pr = identity();
        changeRows(pr, rowsInPs.get(i).rowA, rowsInPs.get(i).rowB);
        double[][] mtilda_i = Ms.get(i - 1);
        changeRows(mtilda_i, rowsInPs.get(i).rowA, rowsInPs.get(i).rowB);
        changeColumns(mtilda_i, rowsInPs.get(i).rowA, rowsInPs.get(i).rowB);
        Mtildas.add(mtilda_i);
        while (i != 1) {
            i--;
            changeColumns(pl, rowsInPs.get(i).rowA, rowsInPs.get(i).rowB); // p_(n-1)*...*p_(k+1)
            changeRows(pr, rowsInPs.get(i).rowA, rowsInPs.get(i).rowB); // p_(k+1) *...*p_(n-1)
            mtilda_i = multiply(multiply(pl, Ms.get(i - 1)), pr); //M~_k = p_(n-1)*...*p_(k+1)* M_k* p_(k+1) *...*p_(n-1)
            Mtildas.add(mtilda_i);
        }
        for (int c = 0; c < Mtildas.size(); c++) { // invM~_1 * invM~_2 *...* invM~_(n-2)
            for (int r = c + 1; r < size; r++) {
                L[r][c] = -Mtildas.get(Mtildas.size() - 1 - c)[r][c];
            }
        }
        L[size - 1][size - 2] = -Ms.get(Ms.size() - 1)[size - 1][size - 2]; // L = invM~_1 * invM~_2 *...* invM~_(n-2)*invM_(n-1)
    }

    static void changeRows(double[][] mat, int r, int k) {
        double[] temp = new double[size];
        for (int i = 0; i < size; i++) {
            temp[i] = mat[k][i];
            mat[k][i] = mat[r][i];
            mat[r][i] = temp[i];
        }
    }

    static void changeColumns(double[][] mat, int c, int k) {
        double[] temp = new double[size];
        for (int i = 0; i < mat.length; i++) {
            temp[i] = mat[i][k];
            mat[i][k] = mat[i][c];
            mat[i][c] = temp[i];
        }
    }

    static void show(double[][] matrix) {
        System.out.println();
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                System.out.print(matrix[i][j] + "\t\t");
            }
            System.out.println();
        }
    }

    static double frobeniusNorm(double mat[][]) {
        double sumSq = 0;
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat[0].length; j++) {
                sumSq += mat[i][j] * mat[i][j];
            }
        }
        return Math.sqrt(sumSq);
    }

    static double[][] identity() { // return I[n][n]
        double[][] a = new double[size][size];
        for (int i = 0; i < size; i++)
            a[i][i] = 1;
        return a;
    }

    static double[][] multiply(double[][] a, double[][] b) { // return A*B
        double[][] c = new double[a.length][b[0].length];
        for (int i = 0; i < a.length; i++)
            for (int j = 0; j < b[0].length; j++)
                for (int k = 0; k < b.length; k++)
                    c[i][j] += a[i][k] * b[k][j];
        return c;
    }

    static double[][] subtract(double[][] a, double[][] b) { // return c = a - b
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                a[i][j] = a[i][j] - b[i][j];
        return a;
    }
}
