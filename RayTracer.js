/* A skeleton of this file was written by Duncan Levear in Spring 2023 for CS3333 at Boston College */


//for the images that render black, the box part of the node is undefined

//am i getting some hit values of t < 0????????????????

//hits of node not calling itself in a loop


import {Vector3, vectorSum, vectorDifference, vectorScaled} from './Vector3.js'

const EPSILON = 1e-9
let recursion_amount = 0


export class RayTracer {
    constructor(sceneInfo, image) {
        this.scene = sceneInfo;
        this.image = image;
        // clear image all white
        for (let i = 0; i < image.data.length; i++) {
            image.data[i] = 255;
        }
    }

    putPixel(row, col, r, g, b) {
        /*
        Update one pixel in the image array. (r,g,b) are 0-255 color values.
        */
        if (Math.round(row) != row) {
            console.error("Cannot put pixel in fractional row");
            return;
        }
        if (Math.round(col) != col) {
            console.error("Cannot put pixel in fractional col");
            return;
        }
        if (row < 0 || row >= this.image.height) {
            return;
        }
        if (col < 0 || col >= this.image.width) {
            return;
        }

        const index = 4 * (this.image.width * row + col);
        this.image.data[index + 0] = Math.round(r);
        this.image.data[index + 1] = Math.round(g);
        this.image.data[index + 2] = Math.round(b);
        this.image.data[index + 3] = 255;
    }

    render() {
        /*
        For every pixel of this.image, compute its color, then use putPixel() to set it. 
        */
        // TODO
        
        this.root = new AABBNode(this.scene.a_geometries)
        
          //set up camera and use to pass to ray
        const [e, u, v, w] = this.setupCamera()
        
        // 1. For Each Pixel
        for (let row=0; row < this.image.height; row++) {
            for (let col=0; col < this.image.width; col++) {
    
                // 2. Compute ray
                const ray = this.pixelToRay(row,col, e,u, v, w)
                
                const color = this.traceRay(ray)
                color.scaleBy(255)
                
                this.putPixel(row,col,color.x, color.y, color.z)
            }
        }
       
       
    }
    

    setupCamera() {
        const e = this.scene.v3_eye
        const eyeOut = this.scene.v3_eyeOut
        const up = this.scene.v3_up
                
        const w = vectorScaled(eyeOut,-1)        
        const u = up.crossProduct(w)
        const v = w.crossProduct(u)
        u.normalize()
        v.normalize()
        w.normalize()
        
        
        return [e, u, v, w]
        
    }
        
    pixelToRay(row, col, e, u, v, w) {
        /*
        Convert from pixel coordinate (row,col) to a viewing ray. Return an instance of the Ray class. 
        */
        // TODO
        
     //   new with camera positioning
        const A = e //camera pt
        const d = this.scene.f_imageplaneDistance
        const w_scaled_d = vectorScaled(w, d)
       
        //compute topleft
        const first_term = vectorDifference(e, w_scaled_d)
        const second = vectorScaled(v, (this.scene.f_imageplaneHeight/2))
        const third = vectorScaled(u, (this.scene.f_imageplaneWidth/2))     //an image rendered when this was v not u
        const first_addition = vectorSum(first_term, second)
        
        const topLeft = vectorDifference(first_addition, third)
        
        //compute ray for every pixel
        const squareHeight = this.scene.f_imageplaneHeight/ this.scene.i_height
        const squareWidth = this.scene.f_imageplaneWidth/ this.scene.i_width
        
        
        const first_pixel_second_term = vectorScaled(v,(.5*squareHeight))
        const first_pixel_third_term = vectorScaled(u,(.5*squareWidth))
        const dif = vectorDifference(topLeft, first_pixel_second_term)
        const firstPixel = vectorSum(dif, first_pixel_third_term)
        
        const B_second_term = vectorScaled(v,(row*squareHeight))
        const B_third_term =  vectorScaled(u,(col*squareWidth))
        
        let B = firstPixel.increaseByMultiple(B_second_term, -1)
        B = B.increaseBy(B_third_term)
        
        const direction = vectorDifference(B,A)
        const ray = new Ray(A,direction)

        return ray
    }
    
    traceRay(ray) {
        /*
        Determine the color of the given ray based on this.scene.
        */
        // TODO
        
        
//        const hits = ray.hitsOfNode(this.root);
        const hits = ray.allHits(this.scene.a_geometries);
        
        let minHit = {t: Infinity}
        for (const h of hits) {
            if (h.t < minHit.t && h.t > EPSILON) {
                minHit = h
            }
        }
        
        if (minHit.t === Infinity) {
            return new Vector3(0,0,0)
        }
        
        const c = this.getColor(minHit)
        return c
        
    }
    
    getColor(record) {
        /*
        Determine the color that results from the given HitRecord.
        */
        // TODO      
        let color_added = new Vector3 (0,0,0)
        let reflectedAmount = 0
        let reflectedColor = new Vector3 (0,0,0)

        
        const lights = this.scene.a_lights
        for (let i = 0; i <lights.length; i++) {
            //(lights[i].v3_position)
            const light = lights[i]
            const diffuse_specular = this.whatLight(record, light)
            color_added.increaseBy(diffuse_specular)
        }
            //check if reflective
            if (record.struckGeometry.j_material.f_reflectance>0 && recursion_amount < 5 ) {
                recursion_amount +=1 
                reflectedAmount = record.struckGeometry.j_material.f_reflectance
                //console.log(reflectedAmount)
                reflectedColor = this.reflected(record)
                recursion_amount -=1
                //console.log("reflected color", reflectedColor)
                
                color_added.increaseBy(reflectedColor.scaleBy(reflectedAmount))

                
            }
    
        return color_added

    }
 
    
    // To add shading, break it into steps: whatLight(), diffuse(), highlight(), or similar
    whatLight(hit, light_source) {
   
        
        //compute shadow and return black if shadow is there SHADOW CODE
        const point = hit.pt
        
        const shadowRay_dir = vectorDifference(light_source.v3_position, point)  
        const shadowRay = new Ray(point, shadowRay_dir)

      //  const shadow_hit = shadowRay.hitsOfNode(this.root)
        
        const shadow_hit = shadowRay.allHits(this.scene.a_geometries)

            
        for (const hit of shadow_hit) { 

            if (hit.t > 0.0001 && hit.t < 1) {

            return new Vector3([0, 0, 0]);
            }
        }

        const light_pos = light_source.v3_position
        const toLight = vectorDifference(light_pos, point).normalize()
        

        const dif_color = this.diffuse(hit, light_source, toLight) 
        
        const spec_light = this.specular(hit, light_source, toLight)
     
        const returnMe = new Vector3(0,0,0)
        returnMe.increaseBy(dif_color)
        returnMe.increaseBy(spec_light)
        return returnMe
    }
    
    
    reflected(hit){
                    
        const mirrorRay = this.bounce(hit)
        
        const mirror_hit = this.traceRay(mirrorRay)
 
        return mirror_hit
        
    }
    
    bounce(hit){
        
        const viewingRay = hit.ray

        const vectored_viewingRay = viewingRay.dir 
        
        const normal = hit.normal
        
        const inverse = new Vector3(vectored_viewingRay).scaleBy(-1)

        const top = inverse.dotProduct(normal)
        const bottom = normal.dotProduct(normal)
        
        const bounced_beginning = new Vector3(normal).scaleBy(2*(top/bottom))
        
        const bounced = vectorDifference(bounced_beginning, inverse)
        
        const bouncedRay = new Ray(hit.pt, bounced)
        //console.log("bounced Ray", bouncedRay)
    
        return bouncedRay
        
    
    }
    

    
    diffuse(hit, light, toLight) { //lambert computation
        
 
        const normal = hit.normal
        const alignment = toLight.dotProduct(normal) 
        //console.log(alignment)
        
        if (alignment<0){
            return new Vector3(0,0,0)
        }
    

        const color = vectorScaled(hit.struckGeometry.j_material.v3_diffuse, alignment)
        
        
        color.scaleBy(light.f_intensity)

       
        return color
        
    }
    
    specular (hit, light_source, toLight) {      //phong computation

        const point = hit.pt

        const specularity_power = hit.struckGeometry.j_material.f_specularity
        
        if (specularity_power<0 || specularity_power===undefined){
            return new Vector3(0,0,0)
        }
        
        const to_eye = vectorDifference(this.scene.v3_eye, point)
        
        //const light_pos = light_source.v3_position
        
        const normal = hit.normal
        
        const alpha = 2*toLight.dotProduct(normal)
        const outgoingLight = vectorDifference(vectorScaled(normal,alpha), toLight)
        outgoingLight.normalize()
        to_eye.normalize()
                
   
       
        let s = outgoingLight.dotProduct(to_eye)
     
        if (s < 0) {
            s = 0
        }
        
        s = Math.pow(s,specularity_power)
                
        const final_spec = new Vector3(1,1,1)
        final_spec.scaleBy(s*light_source.f_intensity)
        
        return final_spec
    }
}

class Ray {
    constructor(start, dir) {
        this.start = start; //A x0,y0,z0
        this.dir = dir; //A + this.dir = B
        
        // b = x1,y1,z1
        // P = start + t*dir
    }


    tToPt(t) {
        const ret = new Vector3(this.start).increaseByMultiple(this.dir, t);
        return ret;
    }
    
    allHits(geometries) {
        let ret = [];
        for (const g of geometries) {
            const record = this.hit(g);
            if (record.length === undefined) {
                console.error("Return type of hit() should be an array.");
            }
            ret = ret.concat(record);
        }
        return ret;
    }
    
    hit(g) {
        if (g.s_type === 'sphere') {
            return this.hitSphere(g);
        }
        else if (g.s_type === 'sheet') {
            return this.hitSheet(g);
        }
        else if (g.s_type === 'box') {
            return this.hitBox(g);
        }
        else if (g.s_type === 'cylinder') {
            return this.hitCylinder(g);
        }
        else if (g.s_type === 'triangle') {
            return this.hitTriangle(g);
        }
        else if (g.s_type === 'cone') {
            return this.hitCone(g);
        }
        else {
            console.error("Shape of type " + g.s_type + " is not supported");
        }
    }
    
    hitSheet(g) {
        /*
        Compute the intersection between the ray (this) and the given geometry g, a sheet.
        Return an instance of the HitRecord class.
        */
        // aliases for shorter typing
        const pt0 = g.v3_pt0;
        const pt1 = g.v3_pt1;
        const pt2 = g.v3_pt2;
        // compute d, normal, edge1, edge2 once only, to save time
        if (g.edge1 === undefined) {
            g.edge1 = vectorDifference(pt0, pt1);
            g.edge2 = vectorDifference(pt2, pt1);

            // edge1 and edge2 assumed to be orthogonal
            const unit1 = vectorDifference(pt0, pt1).normalize();
            const unit2 = vectorDifference(pt2, pt1).normalize();
            if (Math.abs(unit1.dotProduct(unit2)) > 0.01) {
                console.error(`Edges ${edge1} and ${edge2} are not orthogonal`);
            }

            // assume pts listed in ccw order, e.g. [1, 0, 0], [0,0,0], [0, 1, 0]
            g.normal = unit2.crossProduct(unit1);
            g.normal.normalize();

            // ray-plane intersection
            g.d = g.normal.dotProduct(pt1);
        }
        
        const t = (g.d - g.normal.dotProduct(this.start))/g.normal.dotProduct(this.dir);
        const pt = this.tToPt(t);
        
        // check if pt is within sheet
        let alpha = vectorDifference(pt,pt1).dotProduct(g.edge1);
        alpha /= g.edge1.dotProduct(g.edge1);
        let beta = vectorDifference(pt,pt1).dotProduct(g.edge2);
        beta /= g.edge2.dotProduct(g.edge2);

        if (alpha < 0 || alpha > 1 || beta < 0 || beta > 1) {
            // hit doesn't count
            return [];
        }
        
        //maybe return this top one later
        const hit = new HitRecord(this, t, pt, g, g.normal)
        return [hit];
        
        
    }

    hitSphere(g) {
        /*
        Compute the intersection between the ray (this) and the given geometry g, a sphere.
        Return an instance of the HitRecord class.
        */
        // TODO
        const r = g.f_radius        //sphere radius
        
        const u = vectorDifference(this.start, g.v3_center)         //A-q , start - spheres center
        const v = this.dir          //B-A
        const a = v.dotProduct(v)           //<v,v>
        const b = 2* (u.dotProduct(v))
        const c = u.dotProduct(u) - (r*r)
 
        
        const discriminant = (b*b) - (4*a*c)
        if (discriminant < 0) {
            const t = (-b + Math.sqrt(discriminant)) / (2*a)
            const pt = this.tToPt(t)
            const pt_normal = vectorDifference(pt, g.v3_center).normalize()
            
            return []
        }
        

        
        //interesections
        const t1 = (-b + Math.sqrt(discriminant)) / (2*a)
        const t2 = (-b - Math.sqrt(discriminant)) / (2*a)
        
        
        //points where ray hit at t1 and t2
        const pt1 = this.tToPt(t1);
        const pt2 = this.tToPt(t2);
        

        //normals of circle are vector b/w pt - center of circle
        const pt1_normal = vectorDifference(pt1, g.v3_center).normalize()
        const pt2_normal = vectorDifference(pt2, g.v3_center).normalize()

        
        const hit1 = new HitRecord(this, t1, pt1, g, pt1_normal)
        const hit2 = new HitRecord(this, t2, pt2, g, pt2_normal)
      
    
        //return instance of hit record class
        if (t1 < t2) {
            return [hit1, hit2]
        }
        else {
            return [hit2, hit1]
        }
        
        
        }

    hitBox(g) {
        /*
        Compute the intersection between the ray (this) and the given geometry g, a box.
        Return an instance of the HitRecord class.
        */
        // TODO
        
        
        const boxX = [g.v3_minPt.x, g.v3_minPt.x +g.v3_dim.x]
        const boxY = [g.v3_minPt.y, g.v3_minPt.y + g.v3_dim.y]
        const boxZ = [g.v3_minPt.z, g.v3_minPt.z + g.v3_dim.z]
        const e = this.start
        const v = this.dir 
        
        let intervalX, intervalY, intervalZ
        
        //X INTERVAL
        
        //if v.x is 0
        if (v.x == 0) {
            if ((boxX[0])<= e.x <= (boxX[1])) {
                intervalX = [-1* Infinity, Infinity]
            }
            else {
                return []
            }
        }
        
        // if v is not 0
        else {
             intervalX = [(boxX[0] - e.x)/ v.x, (boxX[1] - e.x)/v.x]
            //whap the values if the second is bigger than the second
             if (((boxX[0] - e.x)/ v.x) > ((boxX[1] - e.x)/v.x)) {
                 intervalX = [(boxX[1] - e.x)/v.x, (boxX[0] - e.x)/ v.x]
                }
        }
        
        //Y INTERVAL
        
        //if v.y is 0
        if (v.y == 0) {
            if ((boxY[0])<= e.y <= (boxY[1])) {
                intervalY = [-1* Infinity, Infinity]
            }
            else {
                return []
            }
        }
        
        // if v is not 0
        else {
             intervalY = [(boxY[0] - e.y)/ v.y, (boxY[1] - e.y)/v.y]
            //whap the values if the second is bigger than the second
             if (((boxY[0] - e.y)/ v.y) > ((boxY[1] - e.y)/v.y)) {
                 intervalY = [(boxY[1] - e.y)/v.y, (boxY[0] - e.y)/ v.y]
                }
        }
        
         //Z INTERVAL
        
        //if v.z is 0
        if (v.z == 0) {
            if ((boxZ[0])<= e.z <= (boxZ[1])) {
                intervalZ = [-1* Infinity, Infinity]
            }
            else {
                return []
            }
            
            if (boxZ[0] > boxZ[1]) {
                 intervalZ = [(boxZ[1]), (boxZ[0])]
            }
        }
        
        // if v is not 0
        else {
             intervalZ = [(boxZ[0] - e.z)/ v.z, (boxZ[1] - e.z)/v.z]
            //whap the values if the second is bigger than the second
             if (((boxZ[0] - e.z)/ v.z) > ((boxZ[1] - e.z)/v.z)) {
                 intervalZ = [(boxZ[1] - e.z)/v.z, (boxZ[0] - e.z)/ v.z]
            }
        }
        
        
        //compute intersection of the intervals-- think of using a number line
        const max_start = (Math.max(intervalX[0], intervalY[0], intervalZ[0]))
       
        //have to amke sure the min is after the max
        const min_end = (Math.min(intervalX[1], intervalY[1], intervalZ[1]))
        
        
        //if intervals overlap
        if (max_start <= min_end) {
            const t1 = max_start
            const t2 = min_end

            //points where ray hit at t1 and t2
            const pt1 = this.tToPt(t1);
            const pt2 = this.tToPt(t2);
            
            let normal1
            
            if (Math.abs(pt1.x-g.v3_minPt.x) < EPSILON) {
                normal1 = new Vector3(-1,0,0)
            }
            else if (Math.abs(pt1.x - (g.v3_minPt.x +g.v3_dim.x)) < EPSILON) {
                normal1 = new Vector3(1,0,0)
            }
            else if (Math.abs(pt1.y-g.v3_minPt.y) < EPSILON) {
                normal1 = new Vector3(0,-1,0)
            }
            else if (Math.abs(pt1.y - (g.v3_minPt.y +g.v3_dim.y)) < EPSILON) {
                normal1 = new Vector3(0,1,0)
            }
            if (Math.abs(pt1.z-g.v3_minPt.z) < EPSILON) {
                normal1 = new Vector3(0,0,-1)
            }
            else if (Math.abs(pt1.z - (g.v3_minPt.z +g.v3_dim.z))< EPSILON ) {
                normal1 = new Vector3(0,0,1)
            }
            
            let normal2
            
            if (Math.abs(pt2.x-g.v3_minPt.x) < EPSILON) {
                normal2 = new Vector3(-1,0,0)
            }
            else if (Math.abs(pt2.x - (g.v3_minPt.x +g.v3_dim.x)) < EPSILON) {
                normal2 = new Vector3(1,0,0)
            }
            else if (Math.abs(pt2.y-g.v3_minPt.y) < EPSILON) {
                normal2 = new Vector3(0,-1,0)
            }
            else if (Math.abs(pt2.y - (g.v3_minPt.y +g.v3_dim.y)) < EPSILON) {
                normal2 = new Vector3(0,1,0)
            }
            if (Math.abs(pt2.z-g.v3_minPt.z) < EPSILON) {
                normal2 = new Vector3(0,0,-1)
            }
            else if (Math.abs(pt2.z - (g.v3_minPt.z +g.v3_dim.z))< EPSILON ) {
                normal2 = new Vector3(0,0,1)
            }
     

            const hit1 = new HitRecord(this, t1, pt1, g, normal1)
            const hit2 = new HitRecord(this, t2, pt2, g, normal2)

            return [hit1, hit2]
        }
        else {
            
            return []
        }
    }
    
     getM(a,b,c,d,e,f,g,h,i){
        return (a*(e*i-h*f)) + (b*(g*f-d*i)) + (c*(d*h - e*g))
    }
    
    beta(a,b,c,d,e,f,g,h,i,j,k,l){
        return (j*(e*i - h*f)) + (k*(g*f - d*i)) + (l*(d*h-e*g))
        
    }
    
    gamma(a,b,c,d,e,f,g,h,i,j,k,l){
        return (i*(a*k-j*b)) + (h*(j*c-a*l)) + (g*(b*l-k*c))
    }
    
    getT(a,b,c,d,e,f,g,h,i,j,k,l){
        
        let returnMe = f*(a*k-j*b) + e*(j*c-a*l) + d*(b*l-k*c)
        returnMe = returnMe * -1
        return returnMe
    }
    
    
     hitTriangle(g){
        
        const a = g.v3_pt0;
        const b = g.v3_pt1;
        const c = g.v3_pt2;
        
        const e = this.start
        const v = this.dir 
        
        
        const edge1 = vectorDifference(b,a)
        const edge2 = vectorDifference(b,c)
        
        const norm = edge2.crossProduct(edge1).normalize()
        
        //direction = v
        //e is start 
        
        const m = this.getM(b.x-a.x, b.y-a.y, b.z-a.z, c.x-a.x, c.y-a.y, c.z-a.z, -v.x, -v.y, -v.z)
        
        let beta = this.beta(b.x-a.x, b.y-a.y, b.z-a.z, c.x-a.x, c.y-a.y, c.z-a.z, -v.x, -v.y, -v.z, e.x-a.x, e.y-a.y, e.z-a.z)
        
        beta = beta/m
         
        let gam = this.gamma(b.x-a.x, b.y-a.y, b.z-a.z, c.x-a.x, c.y-a.y, c.z-a.z, -v.x, -v.y, -v.z, e.x-a.x, e.y-a.y, e.z-a.z)
        
        gam = gam/m
         
        let t = this.getT(b.x-a.x, b.y-a.y, b.z-a.z, c.x-a.x, c.y-a.y, c.z-a.z, -v.x, -v.y, -v.z, e.x-a.x, e.y-a.y, e.z-a.z)
        
        t = t/m
        const pt = this.tToPt(t)
        
        const alpha = 1-beta-gam
        
    
        if (alpha<0 || alpha>1){
            return []
        }
        
        if (beta<0 || beta>1){
            return []
        }
         
        if (gam<0 || gam>1){
            return []
        }
         
        const hit =  new HitRecord(this, t, pt, g, norm)
        return [hit]
   
    }
        
        
    hitsOfNode(node) {
        //return list of hits between ray and objects contained at or below the given AABB-tree node
        
        //check if ray hits node's bounding box
        const hits = this.hitBox(node.box)
       
        //if ray doesn't hit node.boundingbox
        if (hits.length == 0) {
            return [] 
        }
        if (node.leaf) {
            return this.allHits(node.geometries)    
        }
        
  
        const leftHits = this.hitsOfNode(node.left)
        
        const rightHits = this.hitsOfNode(node.right)
        
        const fullHits = leftHits.concat(rightHits)
        
        return fullHits

        
        }

    
    hitCylinder(g) {
        const center = g3.v3_center
        const height = g.f_height
        const radius = g.f_radius
        
        const testSphere = {
            v3_center: new Vector3(center.x, 0, center.z),
            f_radius: radius,
        }
        const testRay = new Ray(new Vector3(this.start.x, 0, this.start.z, new Vector3(this.dir.x, 0, this.dir.z)))
        const hits = testRay.hitSphere(testSphere)
        const ret = []
        for (const h of hits) {
            const pt = this.tToPt(h.t)
            if (center.y - height/2 < pt.y && pt.y < center.y + height/2) {
                h.struckGeometry = g
                h.pt = pt
                h.ray = this
                ret.push(h)
            }
        }
        for (const m of [-1,1]) {
            const capCenter = new Vector3(center.x, center.y + m*height/2, center.z)
            const normal = new Vector3(0,m,0)
            const rhs = normal.dotProduct(capCenter)
            const t = (rhs - this.start.dotProduct(normal))/this.dir.dotProduct(normal)
            const hit = new HitRecord(this, t, this.tToPt(t), g, normal)
            const ptToAxis = vectorDifference(capCenter, hit.pt)
            if (ptToAxis.norm() < radius) {
                ret.push(hit)
            }
        }
        return ret
    }
    
    hitCone(g) {
        const px = g.v3_tip.x
        const py = g.v3_tip.y
        const pz = g.v3_tip.z
        const r = g.f_radius
        const h = g.f_height
        
        const ex = this.start.x
        const ey = this.start.y
        const ez = this.start.z
        
        const vx = this.dir.x
        const vy = this.dir.y
        const vz = this.dir.z
        
        const ux = ex - px
        const uy = ey - py
        const uz = ez - pz
        
        
        
        const a = Math.pow(vx,2) + Math.pow(vz,2) - ((Math.pow(r,2)/Math.pow(h,2))*Math.pow(vy,2))
        const b = (2*vx*ux) + (2*vz*uz) -  (2*(Math.pow(r,2)/Math.pow(h,2))*vy*uy)
        const c = Math.pow(ux,2) + Math.pow(uz,2) - ((Math.pow(r,2)/Math.pow(h,2))*Math.pow(uy,2))
        
        
        
        const discriminant = (b*b) - (4*a*c)
  
        const ret = []
        if (discriminant > 0) {
            
            const t1 = (-b + (Math.sqrt(discriminant)))/(2*a)
            const t2 = (-b - (Math.sqrt(discriminant)))/(2*a)
            
            var hit1 = 0
            var hit2 = 0
            
            if (t1> 0) {
                
                const pt1 = this.tToPt(t1);
                const pt1_normal = new Vector3((ex + t1*vx - px),(-(Math.pow(r,2)/Math.pow(h,2))*(ey-(t1*vy)-py)), (ez + t1*vz - pz)).normalize()
                hit1 = new HitRecord(this, t1, pt1, g, pt1_normal)
               
                
            }
            
            if (t2> 0) {
                
                const pt2 = this.tToPt(t2);
                const pt2_normal =  new Vector3((ex + t2*vx - px),(-(Math.pow(r,2)/Math.pow(h,2))*(ey-(t2*vy)-py)), (ez + t2*vz - pz)).normalize()

                hit2 = new HitRecord(this, t2, pt2, g, pt2_normal)
               
                
            }
            
            
             if (t1 < t2) {
                ret.push(hit1, hit2)
            }
            else {
                ret.push(hit2, hit1)
            }
            
        }


        for (const m of [-1,1]) {
            const capCenter = new Vector3(g.v3_tip.x, g.v3_tip.y + h, g.v3_tip.z)
            var normal = new Vector3(0,0,0)
            if (h >0) {
                normal = new Vector3(0,-1,0)
            }
            else {
                normal = new Vector3(0,1,0)
            }
           // normal.normalize()
            
            const rhs = normal.dotProduct(capCenter)
            const t = (rhs - this.start.dotProduct(normal))/this.dir.dotProduct(normal)
            const hit = new HitRecord(this, t, this.tToPt(t), g, normal)
            const ptToAxis = vectorDifference(capCenter, hit.pt)
            if (ptToAxis.norm() < r) {
                ret.push(hit)
            }
        }
        
        return ret
    
    }
    
        
    
}

class HitRecord {
    constructor(ray, t, pt, struckGeometry, normal) {
        this.ray = ray; // ray that was involved
        this.t = t; // t-value of intersection along ray
        this.pt = pt; // vector3, point where the ray hit
        this.struckGeometry = struckGeometry; // object that was hit
        this.normal = normal; // normal vector of struckGeometry at pt
    }
}

class AABBNode {
    constructor(g) {
       
//        const firstNode = this.boxAround(g[0])
//        const secondNode = this.boxAround(g[1])
//        var newBox = this.boxAround2(firstNode, secondNode)
//        if (g.length > 3) {
//            for (var i = 2; i < g.length; i++) {
//                const latestBox = this.boxAround2(this.boxAround(g[i]), newBox)
//                newBox = latestBox
//            }
//        }
        
        let box = null
        for (const geo of g) {
            if (box === null) {
                box = this.boxAround(geo)
            }
            else {
                box = this.boxAround2(this.boxAround(geo), box)
            }
            
        }
    
        
         if (g.length <= 3) {
            //make a leaf
            this.geometries = g
            this.left = null
            this.right = null
            this.box = box
            this.leaf = true
            return
        }
      

        this.box = box
        this.leaf = false
        
        
        
        const boxX = this.box.v3_dim.x
        const boxY = this.box.v3_dim.y
        const boxZ = this.box.v3_dim.z
        
        const maxDim = Math.max(boxX, boxY, boxZ)
        
        var sortedGeometries = []
        
        if (maxDim == boxX) {
            sortedGeometries = g.sort((a,b) => this.boxAround(a).v3_minPt.x - this.boxAround(b).v3_minPt.x);
        }
        else if (maxDim == boxY) {
            sortedGeometries = g.sort((a,b) => this.boxAround(a).v3_minPt.y - this.boxAround(b).v3_minPt.y);
        }
        else {
            sortedGeometries = g.sort((a,b) => this.boxAround(a).v3_minPt.z - this.boxAround(b).v3_minPt.z);
        }
        
     
        const leftHalf = sortedGeometries.slice(0, sortedGeometries.length/2)
        const rightHalf = sortedGeometries.slice(sortedGeometries.length/2, sortedGeometries.length)
        this.left = new AABBNode(leftHalf)
        this.right = new AABBNode(rightHalf)
        
    }
    
    
    //g can only be list of 2 boxes
    boxAround2(firstBox, secondBox) {
        
        //call boxAround() and merge 2 boxes into one big box        
        
        const minX = Math.min(firstBox.v3_minPt.x, secondBox.v3_minPt.x);   //is this minX or just x
        const minY = Math.min(firstBox.v3_minPt.y, secondBox.v3_minPt.y);
        const minZ = Math.min(firstBox.v3_minPt.z, secondBox.v3_minPt.z);
        const maxX = Math.max(firstBox.v3_minPt.x + firstBox.v3_dim.x, secondBox.v3_minPt.x + secondBox.v3_dim.x);
        const maxY = Math.max(firstBox.v3_minPt.y + firstBox.v3_dim.y, secondBox.v3_minPt.y + secondBox.v3_dim.y);
        const maxZ = Math.max(firstBox.v3_minPt.z + firstBox.v3_dim.z, secondBox.v3_minPt.z + secondBox.v3_dim.z);
        
         const ret = {
            s_type: 'box',
            v3_minPt: new Vector3(minX, minY, minZ),
            v3_dim: new Vector3(maxX-minX, maxY-minY, maxZ-minZ),
        };
        
        return ret
        
    }
    
    boxAround(g) {
        /*
        Return box around the given geometry.
        */
        if (g.s_type === 'sphere') {
            return this.boxSphere(g);
        }
        else if (g.s_type === 'sheet') {
            return this.boxSheet(g);
        }
        else if (g.s_type === 'box') {
            return g;
        }
        else if (g.s_type === 'cylinder') {
            return this.boxCylinder(g);
        }
        else if (g.s_type === 'cone') {
            return this.boxCone(g);
        }
        else if (g.s_type === 'triangle') {
            return this.boxTriangle(g);
        }
        else {
            console.error("Shape of type " + g.s_type + " is not supported");
        }
        
    }
    
    //pts is an array of points
   boxPoints(pts) {
        /*
        Helper function: return box containing the given points
        */
        let [minX, minY, minZ] = [Infinity, Infinity, Infinity];
        let [maxX, maxY, maxZ] = [-Infinity, -Infinity, -Infinity];
        for (const p of pts) {
            minX = Math.min(minX, p.x);
            minY = Math.min(minY, p.y);
            minZ = Math.min(minZ, p.z);
            maxX = Math.max(maxX, p.x);
            maxY = Math.max(maxY, p.y);
            maxZ = Math.max(maxZ, p.z);
        }
        const ret = {
            s_type: 'box',
            v3_minPt: new Vector3(minX, minY, minZ),
            v3_dim: new Vector3(maxX-minX, maxY-minY, maxZ-minZ),
        };
        return ret;
    }
    boxSheet(g) {
        /*
        Return box around the given geometry, a sheet.
        */
        return this.boxPoints([g.v3_pt0, g.v3_pt1, g.v3_pt2, vectorSum(g.v3_pt2, vectorDifference(g.v3_pt0, g.v3_pt1))]);
    }
    
    boxCylinder(g) {
        /*
        Return box around the given geometry, a cylinder.
        */
        const ret = {
            s_type: 'box',
            v3_minPt: vectorSum(g.v3_center, new Vector3(-g.f_radius, -g.f_height/2, -g.f_radius)),
            v3_dim: new Vector3(2*g.f_radius, g.f_height, 2*g.f_radius),
        };
        return ret;
        
    }
    
    boxSphere(g) {
    
        const ret = {
            s_type: 'box',
            v3_minPt: new Vector3(g.v3_center.x - g.f_radius, g.v3_center.y - g.f_radius, g.v3_center.z - g.f_radius),
            v3_dim: new Vector3(2*g.f_radius, 2*g.f_radius, 2*g.f_radius)
        };
        return ret;
        
    }
    
    boxTriangle(g) {
        
        /*
        Return box around the given geometry, a triangle.
        */        
        const minX = Math.min(g.v3_pt0.x, g.v3_pt1.x, g.v3_pt2.x);
        const minY = Math.min(g.v3_pt0.y, g.v3_pt1.y, g.v3_pt2.y);
        const minZ = Math.min(g.v3_pt0.z, g.v3_pt1.z, g.v3_pt2.z);
        const maxX = Math.max(g.v3_pt0.x, g.v3_pt1.x, g.v3_pt2.x);
        const maxY = Math.max(g.v3_pt0.y, g.v3_pt1.y, g.v3_pt2.y);
        const maxZ = Math.max(g.v3_pt0.z, g.v3_pt1.z, g.v3_pt2.z);

        const ret = {
            s_type: 'box',
            v3_minPt: new Vector3(minX, minY, minZ),
            v3_dim: new Vector3(maxX-minX, maxY-minY, maxZ-minZ),
        };
        return ret;
        
    }
    
    boxCone(g) {
        
        const ret = {
            s_type: 'box',
            v3_minPt: g.v3_tip,
            v3_dim: new Vector3(2*g.f_radius, g.f_height, 2*g.f_radius),
        };
        return ret;
        
    }
    
    
}
