//TIP To <b>Run</b> code, press <shortcut actionId="Run"/> or
// click the <icon src="AllIcons.Actions.Execute"/> icon in the gutter.
/* Написать класс для работы с матрицами комплексных чисел.
В классе должны быть методы, позволяющие сложить, перемножить, разделить матрицы
(при возможности), транспонировать и вычислить определитель.
К классу должен прилагаться консольный интерфейс, позволяющий создать матрицу и задать операцию.
 */
import java.util.InputMismatchException;
import java.util.Scanner;

class Complex {
    double re, im;

    Complex(double re, double im) {
        this.re = re;
        this.im = im;
    }

    Complex add(Complex other) {
        return new Complex(this.re + other.re, this.im + other.im);
    }

    Complex sub(Complex other) {
        return new Complex(this.re - other.re, this.im - other.im);
    }

    Complex mul(Complex other) {
        return new Complex(
                this.re * other.re - this.im * other.im,
                this.re * other.im + this.im * other.re
        );
    }

    Complex div(Complex other) {
        double denom = other.re * other.re + other.im * other.im;
        if (denom < 1e-12) {
            throw new ArithmeticException("Деление на ноль");
        }
        return new Complex(
                (this.re * other.re + this.im * other.im) / denom,
                (this.im * other.re - this.re * other.im) / denom
        );
    }

    double abs() {
        return Math.sqrt(re * re + im * im);
    }

    boolean isZero() {
        return Math.abs(re) < 1e-12 && Math.abs(im) < 1e-12;
    }

    @Override
    public String toString() {
        if (Math.abs(im) < 1e-12) return String.format("%.2f", re);
        if (Math.abs(re) < 1e-12) return String.format("%.2fi", im);
        return String.format("%.2f%+.2fi", re, im);
    }
}

class Matrix {
    Complex[][] data;
    int rows;
    int cols;

    Matrix(int rows, int cols) {
        if (rows <= 0 || cols <= 0) {
            throw new IllegalArgumentException("Размеры матрицы должны быть положительными");
        }
        this.rows = rows;
        this.cols = cols;
        this.data = new Complex[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                data[i][j] = new Complex(0, 0);
            }
        }
    }

    Matrix findSum(Matrix B) {
        if (this.rows != B.rows || this.cols != B.cols) {
            throw new IllegalArgumentException("Сложение невозможно: размеры " + this.rows + "x" + this.cols +
                            " и " + B.rows + "x" + B.cols + " не совпадают");
        }
        Matrix res = new Matrix(this.rows, this.cols);
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.cols; j++) {
                res.data[i][j] = this.data[i][j].add(B.data[i][j]);
            }
        }
        return res;
    }

    Matrix findComp(Matrix B) {
        if (this.cols != B.rows) {
            throw new IllegalArgumentException("Умножение невозможно");
        }
        Matrix res = new Matrix(this.rows, B.cols);
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < B.cols; j++) {
                Complex sum = new Complex(0, 0);
                for (int k = 0; k < this.cols; k++) {
                    sum = sum.add(this.data[i][k].mul(B.data[k][j]));
                }
                res.data[i][j] = sum;
            }
        }
        return res;
    }

    Matrix findQuotient(Matrix B) {
        if (B.rows != B.cols) {
            throw new IllegalArgumentException("Делитель должен быть квадратной матрицей");
        }
        return this.findComp(B.inverse());
    }

    Matrix transpose() {
        Matrix res = new Matrix(this.cols, this.rows);
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.cols; j++) {
                res.data[j][i] = this.data[i][j];
            }
        }
        return res;
    }

    Matrix inverse() {
        if (this.rows != this.cols) {
            throw new IllegalArgumentException("Обратная матрица только для квадратных");
        }
        int n = this.rows;
        Complex[][] a = new Complex[n][2 * n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                a[i][j] = this.data[i][j];
            }
            for (int j = 0; j < n; j++) {
                a[i][n + j] = new Complex(i == j ? 1 : 0, 0);
            }
        }

        for (int i = 0; i < n; i++) {
            int pivot = i;
            for (int k = i + 1; k < n; k++) {
                if (a[k][i].abs() > a[pivot][i].abs()) pivot = k;
            }
            if (a[pivot][i].isZero()) {
                throw new ArithmeticException("Матрица вырождена, обратной не существует");
            }
            if (pivot != i) {
                Complex[] tmp = a[i]; a[i] = a[pivot]; a[pivot] = tmp;
            }

            Complex diag = a[i][i];
            for (int j = 0; j < 2 * n; j++) {
                a[i][j] = a[i][j].div(diag);
            }

            for (int k = 0; k < n; k++) {
                if (k == i) continue;
                Complex factor = a[k][i];
                if (factor.isZero()) continue;
                for (int j = 0; j < 2 * n; j++) {
                    a[k][j] = a[k][j].sub(factor.mul(a[i][j]));
                }
            }
        }

        Matrix res = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                res.data[i][j] = a[i][n + j];
            }
        }
        return res;
    }

    Complex determinant() {
        if (this.rows != this.cols) {
            throw new IllegalArgumentException("Определитель только для квадратных матриц");
        }
        int n = this.rows;

        Complex[][] a = new Complex[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                a[i][j] = this.data[i][j];
            }
        }

        Complex det = new Complex(1, 0);
        int sign = 1;

        for (int i = 0; i < n; i++) {
            int pivot = i;
            for (int k = i + 1; k < n; k++) {
                if (a[k][i].abs() > a[pivot][i].abs()) pivot = k;
            }
            if (a[pivot][i].isZero()) return new Complex(0, 0);

            if (pivot != i) {
                Complex[] tmp = a[i]; a[i] = a[pivot]; a[pivot] = tmp;
                sign = -sign;
            }

            det = det.mul(a[i][i]);

            for (int k = i + 1; k < n; k++) {
                Complex factor = a[k][i].div(a[i][i]);
                for (int j = i; j < n; j++) {
                    a[k][j] = a[k][j].sub(factor.mul(a[i][j]));
                }
            }
        }
        return sign == 1 ? det : new Complex(-det.re, -det.im);
    }
}

public class Main {
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        boolean running = true;

        while (running) {
            System.out.println("\nДоступные операции:\n" +
                    "1. сложить матрицы\n" +
                    "2. перемножить матрицы\n" +
                    "3. разделить матрицы (A * B^-1)\n" +
                    "4. транспонировать\n" +
                    "5. вычислить определитель\n" +
                    "0. выход");
            System.out.print("Выберите номер операции: ");

            try {
                int command = Integer.parseInt(sc.nextLine().trim());

                switch (command) {
                    case 0:
                        running = false;
                        System.out.println("Выход.");
                        break;

                    case 1: {
                        System.out.print("Введите название первой матрицы: ");
                        String n1 = sc.nextLine();
                        System.out.print("Введите название второй матрицы: ");
                        String n2 = sc.nextLine();

                        Matrix A = readMatrix(sc, n1);
                        Matrix B = readMatrix(sc, n2);

                        System.out.println("Матрица суммы:");
                        printMatrix(A.findSum(B));
                        break;
                    }

                    case 2: {
                        System.out.print("Введите название первой матрицы: ");
                        String n1 = sc.nextLine();
                        System.out.print("Введите название второй матрицы: ");
                        String n2 = sc.nextLine();

                        Matrix A = readMatrix(sc, n1);
                        Matrix B = readMatrix(sc, n2);

                        System.out.println("Матрица произведения:");
                        printMatrix(A.findComp(B));
                        break;
                    }

                    case 3: {
                        System.out.print("Введите название делимого: ");
                        String n1 = sc.nextLine();
                        System.out.print("Введите название делителя: ");
                        String n2 = sc.nextLine();

                        Matrix A = readMatrix(sc, n1);
                        Matrix B = readMatrix(sc, n2);

                        System.out.println("Матрица частного (A * B^-1):");
                        printMatrix(A.findQuotient(B));
                        break;
                    }

                    case 4: {
                        System.out.print("Введите название матрицы: ");
                        String n1 = sc.nextLine();
                        Matrix A = readMatrix(sc, n1);
                        System.out.println("Транспонированная матрица:");
                        printMatrix(A.transpose());
                        break;
                    }

                    case 5: {
                        System.out.print("Введите название матрицы: ");
                        String n1 = sc.nextLine();
                        Matrix A = readMatrix(sc, n1);
                        System.out.println("Определитель: " + A.determinant());
                        break;
                    }

                    default:
                        System.out.println("Неверный номер: " + command);
                }
            } catch (NumberFormatException e) {
                System.out.println("Ошибка: нужно ввести целое число.");
            } catch (InputMismatchException e) {
                sc.nextLine();
                System.out.println("Ошибка: неверный формат ввода.");
            } catch (IllegalArgumentException e) {
                System.out.println("Ошибка: " + e.getMessage());
            } catch (ArithmeticException e) {
                System.out.println("Арифметическая ошибка: " + e.getMessage());
            }
        }

        sc.close();
    }

    static Matrix readMatrix(Scanner sc, String name) {
        System.out.print("Введите количество строк для " + name + ": ");
        int r = sc.nextInt();
        sc.nextLine();

        System.out.print("Введите количество столбцов для " + name + ": ");
        int c = sc.nextInt();
        sc.nextLine();

        Matrix m = new Matrix(r, c);
        System.out.println("Введите элементы матрицы " + name + " (" + r + "x" + c + ").");
        System.out.println("Каждый элемент — два числа: вещественная и мнимая части (например: 3 4 → 3+4i)");

        for (int i = 0; i < m.rows; i++) {
            for (int j = 0; j < m.cols; j++) {
                System.out.print("Элемент [" + i + "][" + j + "]: ");
                double re = readDouble(sc);
                double im = readDouble(sc);
                m.data[i][j] = new Complex(re, im);
            }
        }
        sc.nextLine();
        return m;
    }

    static int readPositiveInt(Scanner sc, String what) {
        while (true) {
            try {
                int v = sc.nextInt();
                sc.nextLine();
                if (v <= 0) {
                    System.out.print("Число " + what + " должно быть > 0, повторите: ");
                    continue;
                }
                return v;
            } catch (InputMismatchException e) {
                sc.nextLine();
                System.out.print("Это не число. Введите " + what + " заново: ");
            }
        }
    }

    static double readDouble(Scanner sc) {
        while (true) {
            try {
                double v = sc.nextDouble();
                return v;
            } catch (InputMismatchException e) {
                sc.nextLine();
                System.out.print("Это не число. Повторите ввод: ");
            }
        }
    }

    static void printMatrix(Matrix m) {
        for (int i = 0; i < m.rows; i++) {
            for (int j = 0; j < m.cols; j++) {
                System.out.printf("%18s ", m.data[i][j]);
            }
            System.out.println();
        }
    }
}