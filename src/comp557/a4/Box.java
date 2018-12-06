package comp557.a4;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * A simple box class. A box is defined by it's lower (@see min) and upper (@see max) corner. 
 */
public class Box extends Intersectable {

	public Point3d max;
	public Point3d min;
	
    /**
     * Default constructor. Creates a 2x2x2 box centered at (0,0,0)
     */
    public Box() {
    	super();
    	this.max = new Point3d( 1, 1, 1 );
    	this.min = new Point3d( -1, -1, -1 );
    }	

	@Override
	public void intersect(Ray ray, IntersectResult result) {
		
		// The min and max intersection points with the box slabs in x, y and z
		double tmin, tmax, tymin, tymax, tzmin, tzmax;
		
		// If we view the box from a negative x direction, max.x becomes min.x, and same thing for y and z
		// So for each viewing direction component we need to make the appropriate check
		if(ray.viewDirection.x>0){
			tmin=(min.x-ray.eyePoint.x)/ray.viewDirection.x;
			tmax=(max.x-ray.eyePoint.x)/ray.viewDirection.x;
		}else{
			tmin=(max.x-ray.eyePoint.x)/ray.viewDirection.x;
			tmax=(min.x-ray.eyePoint.x)/ray.viewDirection.x;
		}
		if(ray.viewDirection.y>0){
			tymin=(min.y-ray.eyePoint.y)/ray.viewDirection.y;
			tymax=(max.y-ray.eyePoint.y)/ray.viewDirection.y;
		}else{
			tymin=(max.y-ray.eyePoint.y)/ray.viewDirection.y;
			tymax=(min.y-ray.eyePoint.y)/ray.viewDirection.y;
		}
		if(ray.viewDirection.z>0){
			tzmin=(min.z-ray.eyePoint.z)/ray.viewDirection.z;
			tzmax=(max.z-ray.eyePoint.z)/ray.viewDirection.z;
		}else{
			tzmin=(max.z-ray.eyePoint.z)/ray.viewDirection.z;
			tzmax=(min.z-ray.eyePoint.z)/ray.viewDirection.z;
		}
		// tmin is the max of the minimums, and tmax is the min of the maximums
		tmin=Double.max(tmin, tymin);
		tmin=Double.max(tmin, tzmin);
		tmax=Double.min(tmax, tymax);
		tmax=Double.min(tmax,  tzmax);
		
		// If there is an intersection, and the intersection is in front of the camera
		if(tmin<tmax && tmin>0){
			// update t, p and material of intersection result
			result.material=material;
			result.t=tmin;
			Point3d p=new Point3d(ray.viewDirection);
			p.scale(tmin);
			p.add(ray.eyePoint);
			result.p=p;
			final double epsilon=1e-9;
			// Check on which face of the box the intersection occurred to get the normal
			if(Math.abs(p.x-min.x)<epsilon){
				result.n=new Vector3d(-1,0,0);
			}else if(Math.abs(p.x-max.x)<epsilon){
				result.n=new Vector3d(1,0,0);
			}else if(Math.abs(p.y-min.y)<epsilon){
				result.n=new Vector3d(0,-1,0);
			}else if(Math.abs(p.y-max.y)<epsilon){
				result.n=new Vector3d(0,1,0);
			}else if(Math.abs(p.z-min.z)<epsilon){
				result.n=new Vector3d(0,0,-1);
			}else if(Math.abs(p.z-max.z)<epsilon){
				result.n=new Vector3d(0,0,1);
			}
		}
	}
}
