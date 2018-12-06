package comp557.a4;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.vecmath.Color3f;
import javax.vecmath.Color4f;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * Simple scene loader based on XML file format.
 */
public class Scene {
    
    /** List of surfaces in the scene */
    public List<Intersectable> surfaceList = new ArrayList<Intersectable>();
	
	/** All scene lights */
	public Map<String,Light> lights = new HashMap<String,Light>();

    /** Contains information about how to render the scene */
    public Render render;
    
    /** The ambient light colour */
    public Color3f ambient = new Color3f();
    
    /** The number of recursive reflection calls */
    public int reflections = 0;

    /** 
     * Default constructor.
     */
    public Scene() {
    	this.render = new Render();
    }
    
    /**
     * renders the scene
     */
    public void render(boolean showPanel) {
 
        Camera cam = render.camera; 
        int w = cam.imageSize.width;
        int h = cam.imageSize.height;
        
        render.init(w, h, showPanel);
        
        for ( int i = 0; i < h && !render.isDone(); i++ ) {
            for ( int j = 0; j < w && !render.isDone(); j++ ) {
            	
            	// List of offsetted rays for antialiasing
            	ArrayList<Ray> rays=new ArrayList<Ray>();
            	
            	// Cast multiple rays distributed in a uniform grid in the current pixel for antialiasing.
            	//The number of columns is the floor of the square root of the number of samples.
            	//The number of rows is the floor of samples divided by columns
            	// And there might be a remainder number of samples.
            	int columns=(int)Math.sqrt(render.samples);
            	int rows=render.samples/columns;
            	int remainder=render.samples-rows*columns;
            	// So for instance, 13 samples will result in 4 rows and 3 columns, but with a remainder of 1.
            	// This remainder will be added to the last row, so it will contain 4 columns instead of 3.
            	for(int k=0;k<rows;k++){//iterate through rows
            		//If we are at the last row, add the remainder samples as columns
            		int numCol=columns;
            		if(k==rows-1)
            			numCol+=remainder;
            		
            		for(int l=0;l<numCol;l++){//iterate through columns
            			// The width of a sub-pixel
                    	double areax=(1.0/(numCol+1));
                    	// The height of a sub-pixel
                    	double areay=(1.0/(rows+1));
                    	// If there is no jitter, the ray is directed at the center of the sub-pixel's area
                    	double jitterx=areax/2;
                    	double jittery=areay/2;
                    	// If the scene has jitter, we direct the ray at a random x and y
                    	// within the sub-pixel's area
                    	if(render.jitter){
                        	jitterx=Math.random()*areax;
                        	jittery=Math.random()*areay;
                    	}
                    	// The ray for sub-pixel (k, l) has x offset k*areax+jitterx-0.5 from the center of the main pixel
                    	// and y offset l*areay+jittery-0.5 from the center of the main pixel
                    	double xOffset=k*areax+jitterx-0.5;
                    	double yOffset=l*areay+jittery-0.5;
                    	// Generate the offsetted ray and add it to the list of rays
                    	Ray ray=new Ray();
                		generateRay(i, j, new double[]{yOffset,xOffset}, cam, ray);
                		rays.add(ray);
            		}
            	}
            	
            	// Contains the average color contribution of the offsetted rays
            	Color4f total=new Color4f();
            	// For each offsetted ray, we calculate the pixel color and add it to the total
            	for(Ray ray:rays){
            		// Calulate the ray's intersection
                	IntersectResult result=new IntersectResult();
                	for(Intersectable surface:surfaceList){
                		surface.intersect(ray, result);
                	}
                	result.n.normalize();//Normalize the intersection's normal
                	//Calculate the pixel color for the ray and add it to the total color
                	Color4f c=calcColorValue(ray, result, reflections);
                	total.add(c);//add to total
            	}
            	
            	// Average the contribution of each offsetted ray (divide by number of rays)
            	total.scale((float) (1.0/render.samples));
            	total.clampMax(1);//No color component can be greater than 1 (255 rgb max)
            	//Convert the color to rgb
            	int r = (int)(255*total.x);
                int g = (int)(255*total.y);
                int b = (int)(255*total.z);
                int a = 255;
                int argb = (a<<24 | r<<16 | g<<8 | b); 
                // update the render image with the average color
                render.setPixel(j, i, argb);
                
            }
        }
        
        // save the final render image
        render.save();
        
        // wait for render viewer to close
        render.waitDone();
        
    }
    
    /**
     * This method calculate's the color at a ray's intersection point. I did not use
     * a static method because I wanted to access some of the scene's parameters.
     * @param ray The intersected ray
     * @param result The result of the intersection
     * @param level The recursive level for reflections
     * @return The color for the ray's intersection
     */
    public Color4f calcColorValue(final Ray ray, final IntersectResult result, int level){
    	//By default, the color will be the background color of the render.
    	Color4f L = new Color4f(render.bgcolor.x, render.bgcolor.y, render.bgcolor.z, 0);
    	//If the ray intersects a surface
    	if(result.t<Double.POSITIVE_INFINITY){
    		// First we add the ambiant light contribution
    		L=new Color4f(ambient.x*result.material.diffuse.x,
    				ambient.y*result.material.diffuse.y,ambient.z*result.material.diffuse.z,0);
    		// Then we sum the diffuse and specular contribution of each light
    		for(String key:lights.keySet()){
    			Light light=lights.get(key);//The current light
    			// We first cast a shadow ray from the intersection point to the current light, and
    			// only continue if the shadow ray does not collide with a surface (no shadow)
        		if(!inShadow(result, light, surfaceList, new IntersectResult(), new Ray())){
        			// l is the normalized vector from the intersection point to the light
        			Vector3d l=new Vector3d(light.from);
        			l.sub(result.p);
        			l.normalize();
        			// Bisector is vector between l and v (v=vector from intersection to ray source)
        			Vector3d bisector=new Vector3d(ray.eyePoint);
        			bisector.sub(result.p);
        			bisector.normalize();
        			bisector.add(l);
        			bisector.normalize();
        			
        			// Now we calculate diffuse light which is k_d*I*max(0, n dot l)
        			// Where I is the light's color scaled by the light's power
        			Color4f diffuse=new Color4f(result.material.diffuse.x*light.color.x,
        					result.material.diffuse.y*light.color.y,result.material.diffuse.z*light.color.z,1);
        			diffuse.scale((float)(light.power*Math.max(0, result.n.dot(l))));
        			L.add(diffuse);// Add diffuse contribution
        			
        			// Now calculate specular lighting with blinn-phong model
        			// k_s*I*max(0, n dot bisector)^p, with I being the light color scaled by its power
        			Color4f specular=new Color4f(result.material.specular.x*light.color.x,
        					result.material.specular.y*light.color.y,result.material.specular.z*light.color.z,1);
        			specular.scale((float)(light.power*Math.pow(Math.max(0, result.n.dot(bisector)),result.material.shinyness)));
        			L.add(specular);// Add specular contribution
        			
        			// For bonus points part 11, I added recursive reflections.
        			// I based myself on page 87 of the textbook
        			if(level>0){// If we are not at the last level of recursion
        				// Cast a reflected ray starting from the intersection point
        				Ray refRay=new Ray();
        				// The direction of the reflected ray is the direction
        				// of the original ray reflected about the normal
        				// obtained by r=d-2(d dot n)*n (from textbook)
        				refRay.viewDirection.set(result.n);
        				refRay.viewDirection.scale(-2*ray.viewDirection.dot(result.n));
        				refRay.viewDirection.add(ray.viewDirection);
        				// The ray's origin is the intersection point plus an epsilon offset
        				// to avoid calculating self-reflections
        				refRay.eyePoint=new Point3d(refRay.viewDirection);
        				refRay.eyePoint.scale(1e-9);
        				refRay.eyePoint.add(result.p);        				// Calculate the intersection of the reflected ray
        				IntersectResult res=new IntersectResult();
        				for(Intersectable surface:surfaceList){
                    		surface.intersect(refRay, res);
                    	}
        				res.n.normalize();//make the normal unit length
        				// If the reflected ray intersects a surface
        				if(res.t<Double.POSITIVE_INFINITY){
        					//Recursively call calculateColor to detect the color at the intersection
        					// of the reflected ray, and scale it by the material's reflection (mirror) coefficient
            				Color4f reflection=calcColorValue(ray, res, level-1);
            				reflection=new Color4f(result.material.mirror.x*reflection.x, result.material.mirror.y*reflection.y,
            						result.material.mirror.z*reflection.z, 1);
            				L.add(reflection);// Add reflection contribution
        				}
        			}
        		}
    		}
    	}
		return L;// Return the contribution of all lights
    }
    
    /**
     * Generate a ray through pixel (i,j).
     * 
     * @param i The pixel row.
     * @param j The pixel column.
     * @param offset The offset from the center of the pixel, in the range [-0.5,+0.5] for each coordinate. 
     * @param cam The camera.
     * @param ray Contains the generated ray.
     */
	public static void generateRay(final int i, final int j, final double[] offset, final Camera cam, Ray ray) {
		// From page 75 of the textbook
		// First we calculate an orthonormal basis (u,v,w) for the camera.
		// w is the opposite of the viewing direction
		Vector3d w=new Vector3d(cam.from);
		w.sub(cam.to);
		// length of w before it is normalized gives the distance to the image plane
		double distance=w.length();
		w.normalize();//normalize w
		// Using the field of view, calculate the top of the image plane
		double top=distance*Math.tan(Math.toRadians(cam.fovy/2.0));
		//ratio between width and height of image plane to get right of plane
		double ratio=(double)cam.imageSize.width/cam.imageSize.height;
		double right=top*ratio;//right=(width/height)*top
		double left=-right;// left is negative of right
		double bottom=-top;// bottom is negative of top
		
		Vector3d v=new Vector3d();
		Vector3d u=new Vector3d();
		// u is the cross product of the up vector and w
		u.cross(cam.up, w);
		u.normalize();//normalize u (in case up is not unit length)
		v.cross(u, w);//v is the cross product of u and w
		
		// Now we calculate the ray's direction
		// -direction=-distance*w+a*u+b*v, where a and b are the position of pixels
		// i and j on the image plane. (See page 75 textbook)
		Vector3d direction=new Vector3d();
		// scale u by a=left+(right-left)(j+0.5+offset)/nx
		// an offset was added to j for antialiasing
		u.scale(left+(right-left)*(j+0.5+offset[1])/cam.imageSize.width);
		// scale v by b=top+(top-bottom)(i+0.5+offset)/ny
		v.scale(bottom+(top-bottom)*(i+0.5+offset[0])/cam.imageSize.height);
		w.scale(distance);
		
		// Now the basis vectors are scaled correctly, so we calculate the ray's direction
		direction.add(u);
		direction.add(v);
		direction.sub(w);
		ray.eyePoint=new Point3d(cam.from);//The ray starts from the camera
		ray.viewDirection=direction;//Set the ray's direction
	}

	/**
	 * Shoot a shadow ray in the scene and get the result.
	 * I modified the provided method to pass a List of intersectables instead of a SceneNode
	 * to prevent a problem in TorusMesh.xml (teacher said it was ok).
	 * @param result Intersection result from raytracing. 
	 * @param light The light to check for visibility.
	 * @param surfaces The list of surfaces
	 * @param shadowResult Contains the result of a shadow ray test.
	 * @param shadowRay Contains the shadow ray used to test for visibility.
	 * 
	 * @return True if a point is in shadow, false otherwise. 
	 */
	public static boolean inShadow(final IntersectResult result, final Light light, final List<Intersectable> surfaces, IntersectResult shadowResult, Ray shadowRay) {
		//Cast a shadow ray from the intersection point to the light
		//The shadow ray's direction is the vector from the intersection to the light
		shadowRay.viewDirection.set(light.from);
		shadowRay.viewDirection.sub(result.p);
		// The shadow ray's origin is the intersection point plus an epsilon offset
		// to avoid detecting self-intersections.
		shadowRay.eyePoint=new Point3d(shadowRay.viewDirection);
		shadowRay.eyePoint.scale(1e-9);
		shadowRay.eyePoint.add(result.p);
		//Check if the shadow ray intersects a surface
		for(Intersectable surface:surfaces)
			surface.intersect(shadowRay, shadowResult);
		// If there is an intersection between the initial intersection point and the light, return true. Otherwise return false
		return shadowResult.t<1 && shadowResult.t>0;
	}    
}
