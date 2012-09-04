import org.jblas.*;
import static org.jblas.DoubleMatrix.*;
import static org.jblas.MatrixFunctions.*;

/**
 * Example of how to implement conjugate gradienst with jblas.
 *
 * Again, the main objective of this code is to show how to use
 * jblas, not how to package/modularize numerical code in Java ;)
 *
 * Closely follows http://en.wikipedia.org/wiki/Conjugate_gradient_method
 */
public class CG {
    public static void main(String[] args) {
	new CG().runExample();
    }

    void runExample() {
	int n = 100;
	double w = 1;
	double lambda = 1e-6;
	DoubleMatrix[] ds = sincDataset(n, 0.1);
	DoubleMatrix A = gaussianKernel(w, ds[0], ds[0]).add( eye(n).mul(lambda) );
	DoubleMatrix x = zeros(n);
	DoubleMatrix b = ds[1];

	cg(A, b, x, lambda);
    }
    

    DoubleMatrix cg(DoubleMatrix A, DoubleMatrix b, DoubleMatrix x, double thresh) {
	int n = x.length;
	DoubleMatrix r = b.sub(A.mmul(x));
	DoubleMatrix p = r.dup();
	double alpha = 0, beta = 0;
	DoubleMatrix r2 = zeros(n), Ap = zeros(n);
	while (true) {
	    A.mmuli(p, Ap);
	    alpha = r.dot(r) / p.dot(Ap);
	    x.addi(p.mul(alpha));
            r.subi(Ap.mul(alpha), r2);
	    double error = r2.norm2();
	    System.out.printf("Residual error = %f\n", error);
	    if (error < thresh)
		break;
	    beta = r2.dot(r2) / r.dot(r);
	    r2.addi(p.mul(beta), p);
	    DoubleMatrix temp = r;
	    r = r2;
	    r2 = temp;
	}
	return x;
    }

    /**
     * Compute the Gaussian kernel for the rows of X and Z, and kernel width w.
     */
    public static DoubleMatrix gaussianKernel(double w, DoubleMatrix X, DoubleMatrix Z) {
	DoubleMatrix d = Geometry.pairwiseSquaredDistances(X.transpose(), Z.transpose());
	return exp(d.div(w).neg());
    }

    /**
     * The sinc function (save version).
     *
     * This version is save, as it replaces zero entries of x by 1.
     * Then, sinc(0) = sin(0) / 1 = 1.
     *
     */
    DoubleMatrix safeSinc(DoubleMatrix x) {
	return sin(x).div(x.add(x.eq(0)));
    }


    /**
     * Create a sinc data set.
     *
     *   X ~ uniformly from -4..4
     *   Y ~ sinc(x) + noise * gaussian noise.
     */
    DoubleMatrix[] sincDataset(int n, double noise) {
	DoubleMatrix X = rand(n).mul(8).sub(4);
	DoubleMatrix Y = safeSinc(X) .add( randn(n).mul(noise) );

	return new DoubleMatrix[] {X, Y};
    }
}