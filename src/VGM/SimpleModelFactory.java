package VGM;



public class SimpleModelFactory {
	
	public static Model createModel(String type,double dist, double sill, double range, double nug){
		Model model=null;

		if ("exponential".equals(type)) {
			
			model= new Exponential(dist, sill, range, nug);
			
        }
        if ("gaussian".equals(type)) {
        	
        	model=new Gaussian (dist, sill, range, nug);
        	
        }
        if ("spherical".equals(type)) {
        	
        	model=new Spherical(dist, sill, range, nug);

        }
        if ("pentaspherical".equals(type)) {
        	
        	model=new Pentaspherical(dist, sill, range, nug);

        }
        if ("linear".equals(type)) {
        	
        	model=new Linear (dist, sill, range, nug);
        }
        if ("circular".equals(type)) {
        	
        	model=new Circolar (dist, sill, range, nug);

        }
        if ("bessel".equals(type)) {
        	
        	model=new Bessel (dist, sill, range, nug);

        }
        if ("periodic".equals(type)) {

        	model=new Periodic (dist, sill, range, nug);
        }
        if ("hole".equals(type)) {
        	
        	model=new Hole (dist, sill, range, nug);

        }
        if ("logaritmic".equals(type)) {
        	
        	model = new Logaritmic(dist, sill, range, nug);

        }
        if ("power".equals(type)) {
        	
        	model = new Power(dist, sill, range, nug);

        }
        if ("spline".equals(type)) {
        	
        	model =new Spline(dist, sill, range, nug);

        }
		
		

		return model;
		
	}

}
