package comp557.a4;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * Class for a plane at y=0.
 * 
 * This surface can have two materials.  If both are defined, a 1x1 tile checker 
 * board pattern should be generated on the plane using the two materials.
 */
public class Plane extends Intersectable {
    
	/** The second material, if non-null is used to produce a checker board pattern. */
	Material material2;
	
	/** The plane normal is the y direction */
	public static final Vector3d n = new Vector3d( 0, 1, 0 );
    
    /**
     * Default constructor
     */
    public Plane() {
    	super();
    }

        
    @Override
    public void intersect( Ray ray, IntersectResult result ) { 
    	// The plane is on y=0, so intersection between ray and and plane will occur
    	// at y=0. So Ey+t*Dy=0, we can isolate t and get t=-Ey/Dy
    	double t=-ray.eyePoint.y/ray.viewDirection.y;
    	
    	// If the intersection is in front of the camera (t>0) and closer than the previous one
    	// then we update the intersection result.
    	if(t>0 && t<result.t){
        	result.n.set(n);//the intersection normal is the normal of the plane
        	result.t=t;// update t.
        	// Calculate the intersection point using t
        	result.p=new Point3d(ray.viewDirection);
        	result.p.scale(t);
        	result.p.add(ray.eyePoint);
        	result.material=material;// Set the material to the plane's first material
        	//If we are not in the +x +y or -x -y quadrants, we use the plane's second material instead
        	int x=(int)Math.floor(result.p.x);
        	int z=(int)Math.floor(result.p.z);
        	if(material2!=null && (x+z)%2!=0)
    			result.material=material2;
    	}
    }
    
}
