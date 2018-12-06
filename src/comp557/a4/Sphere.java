package comp557.a4;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * A simple sphere class.
 */
public class Sphere extends Intersectable {
    
	/** Radius of the sphere. */
	public double radius = 1;
    
	/** Location of the sphere center. */
	public Point3d center = new Point3d( 0, 0, 0 );
    
    /**
     * Default constructor
     */
    public Sphere() {
    	super();
    }
    
    /**
     * Creates a sphere with the request radius and center. 
     * 
     * @param radius
     * @param center
     * @param material
     */
    public Sphere( double radius, Point3d center, Material material ) {
    	super();
    	this.radius = radius;
    	this.center = center;
    	this.material = material;
    }
    
    @Override
    public void intersect( Ray ray, IntersectResult result ) {
    	//I am applying the formula from page 77 of the textbook
    	Vector3d d=new Vector3d(ray.viewDirection);//The ray's direction
    	Vector3d e=new Vector3d(ray.eyePoint);//The ray's origin
    	Vector3d ec=new Vector3d(e);// e-c, where c is the sphere's center
    	ec.sub(center);
    	double dec=d.dot(ec);// the dot product of d and e-c
    	double dd=d.dot(d);// the dot product of d and d
    	double ecec=ec.dot(ec);// the dot product of e-c and e-c
    	// The value inside the square root in the formula
    	double innerRoot=dec*dec-dd*(ecec-radius*radius);
    	// There is only an intersection if the square root is not imaginary, i.e. the inside
    	// of the root is greater or equal than 0. Also we need to check for division by 0 for denominator dd
    	if(innerRoot>=0 && dd!=0){
    		// If there is an intersection, we calculate the two possible values for t
        	double root=Math.sqrt(innerRoot);
        	double t1=(-dec+root)/dd;
        	double t2=(-dec-root)/dd;
        	// We keep the smallest of the two answers
        	if(t2<t1)t1=t2;
        	
        	// If the intersection is in front of the camera (t>0) and closer than the previous one
        	// then we update the intersection result
        	if(t1<result.t&&t1>0){
            	result.t=t1;//update t
            	//Calculate the intersection point of the ray and the sphere using t
            	result.p.set(d);
            	result.p.scale(t1);
            	result.p.add(e);
            	// The normal is the vector going from the center to the intersection point
            	Vector3d normal=new Vector3d(result.p);
            	normal.sub(center);
            	normal.scale(1/radius);//Make the normal unit length
            	result.n=normal;
            	result.material=material;//Put the sphere's material in the intersection result
        	}
    	}
    }
    
}
