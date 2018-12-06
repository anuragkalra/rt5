package comp557.a4;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * 
 */
public class Plane extends Intersectable {
    
	/** The second material, if non-null is used to produce a checker board pattern. */
	Material material2;
	
	/** The plane normal is the y direction */
	public static final Vector3d n = new Vector3d(0, 1, 0);
    
    /**
     * Default constructor
     */
    public Plane() {
    	super();
    }

        
    @Override
    public void intersect(Ray ray, IntersectResult result) {
    
    	Vector3d eye = new Vector3d(ray.eyePoint);
		
		if (ray.viewDirection.dot(n) != 0) {
			double t = -(eye.dot(n))/(ray.viewDirection.dot(n));
			if (t < 0 || t > result.t) return;
			Vector3d tdVEC = new Vector3d(ray.viewDirection);
			
			tdVEC.scale(t);
			Point3d p = new Point3d(eye);
			p.add(tdVEC);
				
			//If a material is defined, choose the color based on where it is in the checker pattern
			Material matInit = null;
			boolean x = (((int) Math.ceil(p.x)) % 2) == 0 ? true : false;
			boolean z = (((int) Math.ceil(p.z)) % 2) == 0 ? true : false;
			if (x == z) { 
				matInit = this.material;
			}
			
			else {
				if (material2 == null) {
					matInit = material;
				}
				else {
					matInit = material2;
				}
			}
				
			result.n.set(n);
			result.p.set(p);
			result.t = t;
			result.material = matInit;
		}
    	
    }
    
}