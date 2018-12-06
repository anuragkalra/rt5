/**
 * 
 */

package comp557.a4;

import java.util.HashMap;
import java.util.Map;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

public class Mesh extends Intersectable {

	/**
	 * Static map storing all meshes by name.
	 */
	public static Map<String,Mesh> meshMap = new HashMap<String,Mesh>();

	/**
	 * Name for this mesh. We can reference this to re-use a polygon soup across meshes.
	 */
	public String name;

	/**
	 * The polygon soup.
	 */
	public PolygonSoup soup;

	public Mesh() {
		super();
		this.name = "";
		this.soup = null;
	}			

	@Override
	public void intersect(Ray ray, IntersectResult result) {	

		double aPoint,bPoint,cPoint,dPoint,ePoint,fPoint,gPoint,hPoint;
		double iPoint,jPoint,kPoint,lPoint,transformM,beta, gamma, t;
		
		for(int[] faceVertex: soup.faceList){		
			Point3d A = soup.vertexList.get( faceVertex[0] ).p;
			Point3d B = soup.vertexList.get( faceVertex[1] ).p;
			Point3d C = soup.vertexList.get( faceVertex[2] ).p;		

			aPoint = A.x - B.x;
			bPoint = A.y - B.y;
			cPoint = A.z - B.z;
			dPoint = A.x - C.x;
			ePoint = A.y - C.y;
			fPoint = A.z - C.z;
			gPoint = ray.viewDirection.x;
			hPoint = ray.viewDirection.y;
			iPoint = ray.viewDirection.z;
			jPoint = A.x - ray.eyePoint.x;
			kPoint = A.y - ray.eyePoint.y;
			lPoint = A.z - ray.eyePoint.z;

			transformM = aPoint * ((ePoint * iPoint) - (hPoint * fPoint)) + bPoint * ((gPoint * fPoint) - (dPoint * iPoint)) + cPoint * ((dPoint * hPoint) - (ePoint * gPoint));
			
			t = - (fPoint*(aPoint*kPoint-jPoint*bPoint) + ePoint*(jPoint*cPoint-aPoint*lPoint) + dPoint*(bPoint*lPoint-kPoint*cPoint)) / transformM;

			boolean checker = t > -0.0001 && t < result.t;
			if(checker){
				gamma = (iPoint*(aPoint*kPoint-jPoint*bPoint) + hPoint*(jPoint*cPoint-aPoint*lPoint) + gPoint*(bPoint*lPoint-kPoint*cPoint)) / transformM;
				
				if(gamma > 0 & gamma < 1){
					beta = (jPoint*(ePoint*iPoint-hPoint*fPoint) + kPoint*(gPoint*fPoint-dPoint*iPoint) + lPoint*(dPoint*hPoint-ePoint*gPoint)) / transformM;
					
					if(beta > 0 && beta < 1-gamma){
						//Intersection point is: p = e + td
						Point3d intersection = new Point3d(ray.viewDirection);
						intersection.scale(t);
						intersection.add(ray.eyePoint); 		
						result.p = intersection;
						result.material = this.material;
						result.t = t;

						//Normal calculation
						Vector3d ba = new Vector3d();
						Vector3d ca = new Vector3d();
						Vector3d n = new Vector3d();
						ba.sub(B,A);
						ca.sub(C,A);
						n.cross(ba, ca);
						n.normalize();
						result.n = n;
					}
				}
			}
		}
	}	
}