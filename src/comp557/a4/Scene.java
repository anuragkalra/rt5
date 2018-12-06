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
 * References:
 * ->Uniform Grid Sampling
 * 	http://www.cs.tut.fi/~dkuva2/Lecture7_SamplingAndInterpolation_v3.pdf
 * 	http://web.cs.wpi.edu/~emmanuel/courses/cs563/S10/talks/wk3_p1_wadii_sampling_techniques.pdf
 * ->Stochastic Sampling
 * 	http://www.cs.cmu.edu/afs/cs/academic/class/15869-f11/www/readings/cook86_sampling.pdf
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
    
    public static double epsilonScaler = 1e-9;

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
            	ArrayList<Ray> rays=new ArrayList<Ray>();
            	
            	int columns = (int) Math.sqrt(render.samples);
            	int rows = render.samples / columns;
            	int remainder = render.samples - (columns * rows);
            	
            	for(int k = 0; k < rows; k++){
            		
            		int numCol = columns;
            		if(k == rows - 1) {
            			numCol += remainder;
            		}
            		
            		for(int l = 0; l < numCol; l++) {//iterate through columns
            			// The width of a sub-pixel
                    	double areaX = (1.0 / (numCol+1));
                    	// The height of a sub-pixel
                    	double areaY = (1.0/ (rows+1));
                    	// If there is no jitter, the ray is directed at the center of the sub-pixel's area
                    	double jitterX = areaX / 2;
                    	double jitterY = areaY / 2;
                    	// If the scene has jitter, we direct the ray at a random x and y
                    	// within the sub-pixel's area
                    	if(render.jitter){
                        	jitterX = Math.random() * areaX;
                        	jitterY = Math.random() * areaY;
                    	}
                    	// The ray for sub-pixel (k, l) has x offset k*areax+jitterx-0.5 from the center of the main pixel
                    	// and y offset l*areay+jittery-0.5 from the center of the main pixel
                    	double xOffset = k * areaX + jitterX - 0.5;
                    	double yOffset = l * areaY + jitterY - 0.5;
                    	// Generate the offsetted ray and add it to the list of rays
                    	Ray ray = new Ray();
                		generateRay(i, j, new double[]{yOffset,xOffset}, cam, ray);
                		rays.add(ray);
            		}
            	}
            	
            	// Contains the average color contribution of the offsetted rays
            	Color4f total = new Color4f();
            	// For each offsetted ray, we calculate the pixel color and add it to the total
            	for(Ray ray: rays){
            		// Calulate the ray's intersection
                	IntersectResult result = new IntersectResult();
                	for(Intersectable surface:surfaceList){
                		surface.intersect(ray, result);
                	}
                	result.n.normalize();//Normalize the intersection's normal
                	//Calculate the pixel color for the ray and add it to the total color
                	Color4f c = calcColorValue(ray, result, reflections);
                	total.add(c);//add to total
            	}
            	
            	// Average the contribution of each offsetted ray (divide by number of rays)
            	total.scale((float) (1.0 / render.samples));
            	total.clampMax(1);//No color component can be greater than 1 (255 rgb max)
            	//Convert the color to rgb
            	int r = (int)(255 * total.x);
                int g = (int)(255 * total.y);
                int b = (int)(255 * total.z);
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
    	Color4f L = new Color4f(render.bgcolor.x, render.bgcolor.y, render.bgcolor.z, 0);
    	
    	if(result.t < Double.POSITIVE_INFINITY){
    		// First we add the ambient light contribution
    		float AMB_X = ambient.x * result.material.diffuse.x;
    		float AMB_Y = ambient.y * result.material.diffuse.y;
    		float AMB_Z = ambient.z * result.material.diffuse.z;
    		L = new Color4f(AMB_X, AMB_Y, AMB_Z, 0);
    		// Then we sum the diffuse and specular contribution of each light
    		for(String key: lights.keySet()){
    			Light currentLight = lights.get(key);

        		if(!inShadow(result, currentLight, surfaceList, new IntersectResult(), new Ray())){

        			Vector3d l = new Vector3d(currentLight.from);
        			l.sub(result.p);
        			l.normalize();

        			Vector3d bisector = new Vector3d(ray.eyePoint);
        			bisector.sub(result.p);
        			bisector.normalize();
        			bisector.add(l);
        			bisector.normalize();
        			

        			Color4f diffuse = new Color4f(result.material.diffuse.x * currentLight.color.x,
        					result.material.diffuse.y * currentLight.color.y, result.material.diffuse.z * currentLight.color.z, 1);
        			diffuse.scale((float)(currentLight.power * Math.max(0, result.n.dot(l))));
        			L.add(diffuse);// Add diffuse contribution
        			

        			Color4f specular=new Color4f(result.material.specular.x*currentLight.color.x,
        					result.material.specular.y * currentLight.color.y,result.material.specular.z * currentLight.color.z,1);
        			specular.scale((float)(currentLight.power * Math.pow(Math.max(0, result.n.dot(bisector)), result.material.shinyness)));
        			L.add(specular);// Add specular contribution
        			
        			if(level > 0){// If we are not at the last level of recursion
        				// Cast a reflected ray starting from the intersection point
        				Ray refRay=new Ray();
        				refRay.viewDirection.set(result.n);
        				refRay.viewDirection.scale(-2*ray.viewDirection.dot(result.n));
        				refRay.viewDirection.add(ray.viewDirection);
        				// The ray's origin is the intersection point plus an epsilon offset
        				// to avoid calculating self-reflections
        				refRay.eyePoint = new Point3d(refRay.viewDirection);
        				refRay.eyePoint.scale(1e-9);
        				refRay.eyePoint.add(result.p);        				// Calculate the intersection of the reflected ray
        				IntersectResult res = new IntersectResult();
        				for(Intersectable surface:surfaceList){
                    		surface.intersect(refRay, res);
                    	}
        				res.n.normalize();//make the normal unit length
        				// If the reflected ray intersects a surface
        				if(res.t < Double.POSITIVE_INFINITY){
        					//Recursively call calculateColor to detect the color at the intersection
        					// of the reflected ray, and scale it by the material's reflection (mirror) coefficient
            				Color4f reflection = calcColorValue(ray, res, level - 1);
            				reflection = new Color4f(result.material.mirror.x * reflection.x, result.material.mirror.y * reflection.y,
            						result.material.mirror.z * reflection.z, 1);
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
		
		Vector3d w = new Vector3d(cam.from);
		w.sub(cam.to);
		// length of w before it is normalized gives the distance to the image plane
		double distance = w.length();
		w.normalize();//normalize w
		// Using the field of view, calculate the top of the image plane
		double top = distance * Math.tan(Math.toRadians(cam.fovy/2.0));
		//ratio between width and height of image plane to get right of plane
		double ratio = (double) cam.imageSize.width / cam.imageSize.height;
		double right = top * ratio;//right=(width/height)*top
		double left = -right;// left is negative of right
		double bottom = -top;// bottom is negative of top
		
		Vector3d v = new Vector3d();
		Vector3d u = new Vector3d();
		// u is the cross product of the up vector and w
		u.cross(cam.up, w);
		u.normalize();//normalize u (in case up is not unit length)
		v.cross(u, w);//v is the cross product of u and w
		
		// Now we calculate the ray's direction
		// -direction=-distance*w+a*u+b*v, where a and b are the position of pixels
		// i and j on the image plane. (See page 75 textbook)
		Vector3d direction = new Vector3d();
		
		u.scale(left + (right - left) * (j + 0.5 + offset[1]) / cam.imageSize.width);
		
		v.scale(bottom + (top - bottom) * (i + 0.5 + offset[0])/cam.imageSize.height);
		w.scale(distance);
		
		direction.add(u);
		direction.add(v);
		direction.sub(w);
		ray.eyePoint=new Point3d(cam.from);
		ray.viewDirection=direction;
	}

	/**
	 * 
	 * @param result Intersection result from raytracing. 
	 * @param light The light to check for visibility.
	 * @param surfaces The list of surfaces we need to calculate all shadows
	 * @param shadowResult Contains the result of a shadow ray test.
	 * @param shadowRay Contains the shadow ray used to test for visibility.
	 * 
	 * @return True if a point is in shadow, false otherwise. 
	 */
	public static boolean inShadow(final IntersectResult result, final Light light, final List<Intersectable> surfaces, IntersectResult shadowResult, Ray shadowRay) {
		shadowRay.viewDirection.set(light.from);
		shadowRay.viewDirection.sub(result.p);
		shadowRay.eyePoint = new Point3d(shadowRay.viewDirection);
		shadowRay.eyePoint.scale(Scene.epsilonScaler);
		shadowRay.eyePoint.add(result.p);
		for(Intersectable surface:surfaces)
			surface.intersect(shadowRay, shadowResult);
		boolean finalComp = shadowResult.t < 1 && shadowResult.t > 0;
		return finalComp;
	}    
}
