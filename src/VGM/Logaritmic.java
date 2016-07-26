package VGM;

public class Logaritmic implements Model{
	
	double [] dist;
	double sill;
	double range;
	double nug;
	
	
	public Logaritmic (double[] dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}
	
	

	@Override
	public double[] result() {
        double[] result = new double[dist.length];
        for (int i = 0; i < dist.length; i++) {
            if (dist[i] != 0.0) {
                result[i] = nug + sill * (Math.log(dist[i] / range));
            }
            //System.out.println(result[i]);
        }
        return result;
	}

}
