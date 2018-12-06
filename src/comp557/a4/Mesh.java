package comp557.a4;

import java.util.HashMap;
import java.util.Map;
import javax.vecmath.Point3d;
import javax.vecmath.Tuple3d;
import javax.vecmath.Vector3d;

public class Mesh extends Intersectable {
	
	/** Static map storing all meshes by name */
	public static Map<String,Mesh> meshMap = new HashMap<String,Mesh>();
	
	/**  Name for this mesh, to allow re-use of a polygon soup across Mesh objects */
	public String name = "";
	
	/**
	 * The polygon soup.
	 */
	public PolygonSoup soup;

	public Mesh() {
		super();
		this.soup = null;
	}			
		
	@Override
	public void intersect(Ray ray, IntersectResult result) {
		// I used the formula from the ray tracingcourse slides
		// Loop through every face of the mesh
		for(int[] indices:soup.faceList){
			// a b and c are the 3 vertices of the current face.
			Vector3d a=new Vector3d(soup.vertexList.get(indices[0]).p);
			Vector3d b=new Vector3d(soup.vertexList.get(indices[1]).p);
			Vector3d c=new Vector3d(soup.vertexList.get(indices[2]).p);
			Vector3d bma=new Vector3d(b);bma.sub(a);//b-a
			Vector3d cmb=new Vector3d(c);cmb.sub(b);//c-b
			Vector3d amc=new Vector3d(a);amc.sub(c);//a-c
			
			// The normal is (b-a)cross(c-b)
			Vector3d normal=new Vector3d();
			normal.cross(bma, cmb);
			
			// p is the ray's origin, d is the viewing direction
			Vector3d amp=new Vector3d(a);//a-p
			amp.sub(ray.eyePoint);
			// t = ((a-p) dot n)/(n dot d)
			double t=amp.dot(normal)/(normal.dot(ray.viewDirection));
			// x is the intersection point of the ray on the triangle's plane
			Point3d x=new Point3d(ray.viewDirection);
			x.scale(t);
			x.add(ray.eyePoint);
			
			Vector3d xma=new Vector3d(x);xma.sub(a);//x-a
			Vector3d xmb=new Vector3d(x);xmb.sub(b);//x-b
			Vector3d xmc=new Vector3d(x);xmc.sub(c);//x-c
			
			xma.cross(bma, xma);//(b-a)cross(x-a)
			xmb.cross(cmb,  xmb);//(c-b)cross(x-b)
			xmc.cross(amc, xmc);// (a-c)cross(x-c)
			
			// if (b-a)cross(x-a) >0 and (c-b)cross(x-b)>0 and (a-c)cross(x-c)>0
			// Then the ray's intersection is within the triangle on the plane.
			// If this is the case, and the intersection is in front of camera and closer than before
			// then we update the intersection result.
			if(xma.dot(normal)>0 && xmb.dot(normal)>0 && xmc.dot(normal)>0 && t>0 && t<result.t){
				//set t, the normal, the point and the material of the intersection result
				result.n=normal;
				result.p=x;
				result.material=material;
				result.t=t;
			}
			
		}
		
		
	}
	public static double area(Tuple3d a, Tuple3d b, Tuple3d c){
		return 0.5*Math.abs((a.x-c.x)*(b.y-a.y)-(a.x-b.x)*(c.y-a.y));
	}
}
