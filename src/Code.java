public class Code {


    private int K;
    private double a;
    private double n;
    private double hMin;
    private double hMax;
    private double E;
    private double t1;
    private double r;
    private double x[];
    private double y[];
    private double X[];
    private double Y[];
    private double t[];

    Code() {
        x = new double[1000];
        y = new double[1000];
        t = new double[1000];
        this.K = 50;
        this.x[0] = 1.;
        this.y[0] = 0.;
        this.t[0] = 0.;
        this.a = 0.1;
        this.t1 = 5./1001.;
        this.hMin = this.t1/this.K;
        this.r = this.t1/3;
        this.hMax = this.r*this.hMin;
        this.E = 0.0001;
    }
    public double myFuncX(double t,double x, double y) {
        return -501*x + 500*y;
    }
    public double myFuncY(double t, double x, double y) {
        return 500*x - 501*y;
    }

    public double exactFuncX(double t,double x,double y) {
        return 0.5*(x-y)*Math.pow(Math.E,-1001*t) + 0.5*(x+y)*Math.pow(Math.E,-t);
    }

    public double exactFuncY(double t,double x,double y) {
        return -0.5*(x-y)*Math.pow(Math.E,-1001*t) + 0.5*(x+y)*Math.pow(Math.E,-t);
    }
//    private void MPI(int i) {
//        double x;
//        double y;
//        do {
//            if (Math.abs((x - this.X[i]) + Math.abs(y - this.Y[i])) < this.E) {
//                this.X[i] = x;
//                this.Y[i] = y;
//                break;
//            }else{
//                this.X[i] = x;
//                this.Y[i] = y;
//            }
//        }while (true);
//    }
    private void rungeKutte() {

        double k1, k2, k3;
        this.clear();
        for (int i = 1; i <= 3; i++) {

            k1 = this.myFuncX(this.t[i - 1], this.x[i - 1], this.y[i - 1]);
            k2 = this.myFuncX(this.t[i - 1] + (this.hMin / 3.),
                                 this.x[i - 1],
                              this.y[i - 1] + (this.hMin / 3.) * k1);
            k3 = this.myFuncX(this.t[i - 1] + ((2. / 3.) * this.hMin),
                                 this.x[i - 1],
                              this.y[i - 1] + ((2. / 3.) * this.hMin * k2));
            this.x[i] = this.x[i - 1] + (this.hMin / 4.) * (k1 + (3 * k3));

            k1 = this.myFuncY(this.t[i - 1], this.x[i - 1], this.y[i - 1]);
            k2 = this.myFuncY(this.t[i - 1] + (this.hMin / 3.),
                                 this.x[i - 1],
                              this.y[i - 1] + (this.hMin / 3.) * k1);
            k3 = this.myFuncY(this.t[i - 1] + ((2 / 3.) * this.hMin),
                                 this.x[i - 1],
                              this.y[i - 1] + ((2 / 3.) * this.hMin * k2));
            this.y[i] = this.y[i - 1] + (this.hMin / 4.) * (k1 + (3 * k3));
            System.out.println("x= " +this.x[i] + " y= " + this.y[i] +"   exact X= " + this.exactFuncX(this.t[i],this.x[i],this.y[i]) + " exact Y= " + this.exactFuncY(this.t[i],this.x[i],this.y[i]) + "\n");
            this.t[i] =this.a + i*this.hMin;
        }
    }
    private void clear() {
        clearMass(this.x);
        clearMass(this.y);
        clearMass(this.t);
        this.x[0] = 1.;
        this.y[0] = 0.;
        this.t[0] = 0.1;
    }

    private void clearMass(double[] mass) {
        for (int i = 1; i < mass.length; i++) {
            mass[i] = 0;
        }
    }

    public void process() {
        double x,y;
        this.rungeKutte();
        System.out.println(this.t1);
        System.out.println(this.hMin);
        for(int i = 4;this.t[i-1]<=this.t1;i++){
            this.x[i] = this.x[i - 1] + (this.hMin / 24.) * (9 * this.myFuncX(this.t[i - 1],
                    this.x[i - 1],
                    this.y[i - 1])
                    + (19 * this.myFuncX(this.t[i - 2],
                    this.x[i - 2],
                    this.y[i - 2])) -
                    (5 * this.myFuncX(this.t[i - 3],
                            this.x[i - 3],
                            this.y[i - 3])) +
                    this.myFuncX(this.t[i - 4],
                            this.x[i - 4],
                            this.y[i - 4]));
            this.y[i] = this.y[i - 1] + (this.hMin / 24.) * (9 * this.myFuncY(this.t[i - 1],
                    this.x[i - 1],
                    this.y[i - 1])
                    + (19 * this.myFuncY(this.t[i - 2],
                    this.x[i - 2],
                    this.y[i - 2])) -
                    (5 * this.myFuncY(this.t[i - 3],
                            this.x[i - 3],
                            this.y[i - 3])) + this.myFuncY(this.t[i - 4],
                    this.x[i - 4],
                    this.y[i - 4]));
            System.out.println("x= " +this.x[i] + " y= " + this.y[i] +"   exact X= " + this.exactFuncX(this.t[i],this.x[i],this.y[i]) + " exact Y= " + this.exactFuncY(this.t[i],this.x[i],this.y[i]) + "t= " + this.t[i-1]+"\n");
            this.t[i] =this.a + i*this.hMin;
        }
    }

    public static void main(String[] args) {
        Code obj = new Code();
        obj.process();
    }
}