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
    
    /** The ambient light color */
    public Color3f ambient = new Color3f();
    
    /** Number of reflections for recursion parameter */
    public int reflections = 0;
    
    public static double epsilonScaler = 1e-9;

    public static float totalAccum = 0;
    
    /** 
     * Default constructor.
     */
    public Scene() {
    	this.render = new Render();
    }
    
    public void colorAccumHelper() {
    	float ambientX = 0;
    	float ambientY = 0;
    	float ambientZ = 0;
    	for(int i = 0; i < 5; i++) {
    		ambientX += 0.1;
    		ambientY += 0.1;
    		ambientZ += 0.1;
    	}
    	totalAccum += ambientX + ambientY + ambientZ;
    }
    
    
    /**
     * renders the scene
     */
    public void render(boolean showPanel) {
    	
        Camera cam = render.camera; 
        int w = cam.imageSize.width;
        int h = cam.imageSize.height;
        colorAccumHelper();
        render.init(w, h, showPanel);
        
        for ( int i = 0; i < h && !render.isDone(); i++ ) {
            for ( int j = 0; j < w && !render.isDone(); j++ ) {
            	ArrayList<Ray> allRays = new ArrayList<Ray>();
            	
            	int columns = (int) Math.sqrt(render.samples);
            	int rows = render.samples / columns;
            	int remainder = render.samples - (columns * rows);
            	
            	for(int k = 0; k < rows; k++){
            		
            		int numCol = columns;
            		if(k == rows - 1) {
            			numCol += remainder;
            		}
            		
            		for(int l = 0; l < numCol; l++) {
                    	double areaX = (1.0 / (numCol+1));
                    	double areaY = (1.0/ (rows+1));
                    	double jitterX = areaX / 2;
                    	double jitterY = areaY / 2;
                    	
                    	if(render.jitter == true){
                        	jitterY = Math.random() * areaY;
                        	jitterX = Math.random() * areaX;
                    	}
                    	
                    	double yOffset = l * areaY + jitterY - 0.5;
                    	double xOffset = k * areaX + jitterX - 0.5;

                    	Ray ray = new Ray();
                    	double[] offsetTemp = {yOffset, xOffset};
                		generateRay(i, j, offsetTemp, cam, ray);
                		allRays.add(ray);
            		}
            	}
            	
            	
            	Color4f totalContribution = new Color4f();
            	
            	for(Ray ray: allRays){
            		
                	IntersectResult result = new IntersectResult();
                	for(Intersectable surface:surfaceList){
                		surface.intersect(ray, result);
                	}
                	result.n.normalize();
                	
                	Color4f c = calcColorValue(ray, result, reflections);
                	totalContribution.add(c);
            	}
            	
            	
            	totalContribution.scale((float) (1.0 / render.samples));
            	//need to clamp because of cosine rule
            	totalContribution.clampMax(1);
            	
            	int r = (int)(255 * totalContribution.x);
                int g = (int)(255 * totalContribution.y);
                int b = (int)(255 * totalContribution.z);
                int a = 255;
                int argb = (a<<24 | r<<16 | g<<8 | b); 

                render.setPixel(j, i, argb);
                
            }
        }
        
        // save the final render image
        render.save();
        
        // wait for render viewer to close
        render.waitDone();
        
    }
    
    
    
    //calculates total color contribution value using slide 22 in slides
    public Color4f calcColorValue(final Ray ray, final IntersectResult result, int recLevel){
    	Color4f L = new Color4f(render.bgcolor.x, render.bgcolor.y, render.bgcolor.z, 0);
    	
    	if(result.t < Double.POSITIVE_INFINITY){
    		float AMB_X = ambient.x * result.material.diffuse.x;
    		float AMB_Y = ambient.y * result.material.diffuse.y;
    		float AMB_Z = ambient.z * result.material.diffuse.z;
    		L = new Color4f(AMB_X, AMB_Y, AMB_Z, 0);
    		
    		for(String key: lights.keySet()){
    			Light currentLight = lights.get(key);

        		if(inShadow(result, currentLight, surfaceList, new IntersectResult(), new Ray()) == false){

        			Vector3d l = new Vector3d(currentLight.from);
        			l.sub(result.p);
        			l.normalize();

        			Vector3d bisector = new Vector3d(ray.eyePoint);
        			
        			bisector.sub(result.p);
        			bisector.normalize();
        			bisector.add(l);
        			bisector.normalize();
        			
        			float DIFF_X = result.material.diffuse.x * currentLight.color.x;
        			float DIFF_Y = result.material.diffuse.y * currentLight.color.y;
        			float DIFF_Z = result.material.diffuse.z * currentLight.color.z;
        			
        			Color4f diffuseComponent_i = new Color4f(DIFF_X, DIFF_Y, DIFF_Z, 1);

        			diffuseComponent_i.scale((float)(currentLight.power * Math.max(0, result.n.dot(l))));
        			L.add(diffuseComponent_i);
        			
        			float SPEC_X = result.material.specular.x*currentLight.color.x;
        			float SPEC_Y = result.material.specular.y * currentLight.color.y;
        			float SPEC_Z = result.material.specular.z * currentLight.color.z;
        			
        			Color4f specularComponent_i = new Color4f(SPEC_X, SPEC_Y, SPEC_Z, 1);
        			
        			specularComponent_i.scale((float)(currentLight.power * Math.pow(Math.max(0, result.n.dot(bisector)), result.material.shinyness)));
        			L.add(specularComponent_i);
        			
        			if(recLevel > 0){
        				Ray reflectedRay = new Ray();
        				reflectedRay.viewDirection.set(result.n);
        				reflectedRay.viewDirection.scale(-2 * ray.viewDirection.dot(result.n));
        				reflectedRay.viewDirection.add(ray.viewDirection);
        				
        				reflectedRay.eyePoint = new Point3d(reflectedRay.viewDirection);
        				reflectedRay.eyePoint.scale(Scene.epsilonScaler);
        				reflectedRay.eyePoint.add(result.p);        				
        				IntersectResult res = new IntersectResult();
        				for(Intersectable surface:surfaceList){
                    		surface.intersect(reflectedRay, res);
                    	}
        				
        				res.n.normalize();
        				
        				if(res.t < Double.POSITIVE_INFINITY){
            				Color4f reflection = calcColorValue(ray, res, recLevel - 1);
            				float REF_X_ACCUM = result.material.mirror.x * reflection.x;
            				float REF_Y_ACCUM = result.material.mirror.y * reflection.y;
            				float REF_Z_ACCUM = result.material.mirror.z * reflection.z;
            				reflection = new Color4f(REF_X_ACCUM, REF_Y_ACCUM, REF_Z_ACCUM, 1);            				
            				L.add(reflection);
        				}
        			}
        		}
    		}
    	}
		return L;
    }
    
    
	public static void generateRay(final int i, final int j, final double[] offset, final Camera cam, Ray ray) {
		//Referenced CGPP page 77 for consultation on implementation
		ray.eyePoint = new Point3d(cam.from);
		Vector3d w = new Vector3d(cam.from);
		w.sub(cam.to);
		
		double d = w.length();
		w.normalize();//normalize w
		
		double top = d * Math.tan(Math.toRadians(cam.fovy/2.0));
		
		double aspectRatio = (double) cam.imageSize.width / cam.imageSize.height;
		double right = top * aspectRatio;
		double left = -right;
		double bottom = -top;
		
		Vector3d v = new Vector3d();
		Vector3d u = new Vector3d();
		
		u.cross(cam.up, w);
		u.normalize();
		v.cross(u, w);
		
		
		Vector3d direction = new Vector3d();
		
		u.scale(left + (right - left) * (j + 0.5 + offset[1]) / cam.imageSize.width);
		
		v.scale(bottom + (top - bottom) * (i + 0.5 + offset[0])/cam.imageSize.height);
		w.scale(d);
		
		direction.add(u);
		direction.add(v);
		direction.sub(w);

		ray.viewDirection = direction;
	}

	public static boolean inShadow(final IntersectResult result, final Light light, final List<Intersectable> surfaces, IntersectResult shadowResult, Ray shadowRay) {
		shadowRay.viewDirection.set(light.from);
		shadowRay.viewDirection.sub(result.p);
		
		shadowRay.eyePoint = new Point3d(shadowRay.viewDirection);
		shadowRay.eyePoint.scale(Scene.epsilonScaler);
		shadowRay.eyePoint.add(result.p);
		
		for(Intersectable surface:surfaces) {
			surface.intersect(shadowRay, shadowResult);
		}
		
		boolean finalComp = (shadowResult.t < 1 && shadowResult.t > 0);
		return finalComp;
	}    
}
