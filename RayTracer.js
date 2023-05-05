import {Vector3, vectorSum, vectorDifference, vectorScaled} from './Vector3.js'

const EPSILON = 0.000001;
let curDepth = 0;
export class RayTracer {
    constructor(sceneInfo, image) {
        this.scene = sceneInfo;
        this.image = image;
        // clear image all white
        for (let i = 0; i < image.data.length; i++) {
            image.data[i] = 255;
        }
        this.setCamera();
        this.root = new AABBNode(this.scene.a_geometries);
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
    
    setCamera(){
        this.eyeOutReversed = new Vector3(this.scene.v3_eyeOut).normalize().scaleBy(-1);
        this.eyeRight = this.scene.v3_up.crossProduct(this.eyeOutReversed);
        this.eyeRight.normalize();
        this.eyeUp = this.eyeOutReversed.crossProduct(this.eyeRight);
        this.eyeUp.normalize();
    }

   
    render() {
         /*
        For every pixel of this.image, compute its color, then use putPixel() to set it. 
        */
        // TODO
        for (let row=0; row < this.image.height; row++) {
            for (let col=0; col < this.image.width; col++) {
                this.debug = (
                row == Math.round((0.45)*this.image.height) 
                && col == Math.round((0.48)*this.image.width)
                );
                const ray = this.pixelToRay(row,col);
                const c = this.traceRay(ray);
                c.scaleBy(255);
                c.clampAll(0,255);
                c.roundAll();
                curDepth -= 1;
                this.putPixel(row, col, c.x, c.y, c.z);
            }
        }

    }
    
    pixelToRay(row, col) {
        /*
        use camera params to compute ray direction for a given pixel
        Convert from pixel coordinate (row,col) to a viewing ray. Return an instance of the Ray class. 
        */
        // TODO
    
        const A = this.scene.v3_eye;
        
        const topLeft = new Vector3(this.scene.v3_eye);
        topLeft.increaseByMultiple(this.eyeOutReversed, -this.scene.f_imageplaneDistance);
        topLeft.increaseByMultiple(this.eyeUp, this.scene.f_imageplaneHeight/2);
        topLeft.increaseByMultiple(this.eyeRight, -this.scene.f_imageplaneWidth/2);
        
        const squareHeight = this.scene.f_imageplaneHeight/this.scene.i_height;
        const squareWidth = this.scene.f_imageplaneWidth/this.scene.i_width;
        
        const f = new Vector3(topLeft);
        f.increaseByMultiple(this.eyeUp, -1/2*squareHeight);
        f.increaseByMultiple(this.eyeRight, 1/2*squareWidth);
        
        const B = new Vector3(f);
        B.increaseByMultiple(this.eyeUp, -row*squareHeight);
        B.increaseByMultiple(this.eyeRight, col*squareWidth);
        
        const dir = vectorDifference(B,A);
        
        return new Ray(A, dir);
    }
    
    traceRay(ray) {
        /*
        Determine the color of the given ray based on this.scene.
        
        Compute color of a given ray by computing first hit and calling getColor()
        */
        // TODO
        
        curDepth += 1;
        const hits = ray.hitsOfNode(this.root);

        let T = {t: Infinity};
        for (const h of hits){
            if (h.t < T.t && h.t > EPSILON){
                T = h;
            }      
        }  
        if (T.t === Infinity){
            return new Vector3(0,0,0);
        }
        const color = this.getColor(T);
        return color;
    }
    
    getColor(record) {
        /*
        Determine the color that results from the given HitRecord. Calls WhatLight() 
        */
        // TODO
        let reflection = new Vector3(0,0,0);
        let reflectance = record.struckGeometry.j_material.f_reflectance;
        
        if (reflectance !== undefined && curDepth < 7) {
            reflection = this.reflected(record); 
        }
        else {
            reflectance = 0;
        }

        const ret = new Vector3(0,0,0);
        for (const l of this.scene.a_lights){
            const light = this.whatLight(record,l);
            ret.increaseBy(light);
        }
        return ret.increaseBy(vectorScaled(reflection, reflectance));
    }
    
    bounce(record) {
        
        const start = new Vector3(record.pt)
        const n = new Vector3(record.normal);
        const t = vectorScaled((new Vector3(record.ray.dir)), -1);

        const alpha = 2 * t.dotProduct(n) / n.dotProduct(n);
        const mirrorDir = vectorDifference(vectorScaled(n, alpha), t);
        
        return new Ray(start, mirrorDir);
        
    }
    
    reflected(record){
        const mirrorRay = this.bounce(record);
        const reflection = this.traceRay(mirrorRay);
        curDepth -= 1;
        return reflection; 
    }
    
    whatLight(record,light){
        // compute color contributed by a given light for a given hitRecord
        // code for a single light 
         
        // check for light eligibility (shadows)
        const toLight = vectorDifference(light.v3_position,record.pt);
        const shadowRay = new Ray(record.pt, toLight);
        const hits = shadowRay.hitsOfNode(this.root);
        for (const h of hits) {
            if (EPSILON < h.t && h.t < 1){
                return new Vector3(0,0,0);
            }
        }
        const diffuseColor = this.diffuse(record, light);
        const highlight = this.highlight(record,light);
        let ret = vectorSum(diffuseColor,highlight);
        
        return ret;
    }
    
    diffuse(record, light){
        // compute diffuse component of color from given light. This is the lambert shading
        const toLight = vectorDifference(light.v3_position,record.pt).normalize();
        const weight = toLight.dotProduct(record.normal);
        if (weight < 0){
            return new Vector3(0,0,0);
        } 
        const ret = vectorScaled(record.struckGeometry.j_material.v3_diffuse, weight);
        ret.scaleBy(light.f_intensity);
        return ret;
    }
    
    highlight(record, light){
        // phong shading
        const toLight = vectorDifference(light.v3_position,record.pt).normalize();
        const n = record.struckGeometry.j_material.f_specularity;
        if (n === undefined || n === -1){
            return new Vector3(0,0,0);
        }
        
        const pt = new Vector3(record.pt);
        
        const eye = new Vector3(this.scene.v3_eye);
        
        let toEye = vectorDifference(eye,pt);
        
        const normal = new Vector3(record.normal);
        
        let alpha = 2 * normal.dotProduct(toLight);
        alpha /= normal.dotProduct(normal);
    
        let outgoingLight = vectorDifference(vectorScaled(normal,alpha),toLight);
        toEye = toEye.normalize();
        outgoingLight = outgoingLight.normalize();
                               
        let specularAlignment = toEye.dotProduct(outgoingLight);
        
        if (specularAlignment < 0) {
            specularAlignment = 0;
        }
        
        let s = Math.pow(specularAlignment, n);
        s *= light.f_intensity;
        
        const contribution = new Vector3(1,1,1);
        contribution.scaleBy(s);
        return contribution;
    }
}

class Ray {
    constructor(start, dir) {
        this.start = start;
        this.dir = dir;
    }

    tToPt(t) {
        const ret = new Vector3(this.start).increaseByMultiple(this.dir, t);
        return ret;
    }
    
    allHits(geometries) {
        /* compute all hits between a ray and a list of geometries. Return a list of HitRecords*/
        
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
        else {
            console.error("Shape of type " + g.s_type + " is not supported");
        }
    }
    
    hitSheet(g) {
        /*
        Compute the intersection between the ray (this) and the given geometry g, a sheet.
        Return an instance of the HitRecord class.
        */
    
        const pt0 = g.v3_pt0;
        const pt1 = g.v3_pt1;
        const pt2 = g.v3_pt2;

        if (g.edge1 === undefined) {
            g.edge1 = vectorDifference(pt0, pt1);
            g.edge2 = vectorDifference(pt2, pt1);

            // edge1 and edge2 assumed to be orthogonal
            const unit1 = vectorDifference(pt0, pt1).normalize();
            const unit2 = vectorDifference(pt2, pt1).normalize();
            if (Math.abs(unit1.dotProduct(unit2)) > 0.01) {
                console.error(`Edges ${edge1} and ${edge2} are not orthogonal`);
            }

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
        return [new HitRecord(this, t, pt, g, g.normal)];
    }

    hitSphere(g) {
        /*
        Compute the intersection between the ray (this) and the given geometry g, a sphere.
        Return an instance of the HitRecord class.
        */
       
        // compute u, v, a, b, c
        const center = g.v3_center;
        const start = this.start;
        const r = g.f_radius;
        const u = vectorDifference(start,center)
        const v = this.dir;
        const a = v.dotProduct(v);
        const b = v.dotProduct(u)*2;
        const c = u.dotProduct(u)-r*r;
        
        const x = b*b-4*a*c;
        
        if (x < 0) {
            return [];
        }
        
        let t0 = (-b-Math.sqrt(b*b-4*a*c))/(2*a);
        let t1 = (-b+Math.sqrt(b*b-4*a*c))/(2*a);
         
        const ret = [];
        for (const t of [t0,t1]){
            const pt = this.tToPt(t);
            const normal = vectorDifference(pt, center).normalize();
            ret.push(new HitRecord(this,t,pt,g,normal));
        }
        return ret;
    }

    hitBox(g) {
        /*
        Compute the intersection between the ray (this) and the given geometry g, a box.
        Return an instance of the HitRecord class.
        */
        const e = this.start;
        const v = this.dir;
        const width = g.v3_dim.x;
        const height = g.v3_dim.y;
        const depth = g.v3_dim.z;
        const minPt = g.v3_minPt;
        const boxX = [minPt.x,minPt.x+width];
        const boxY = [minPt.y,minPt.y+height];
        const boxZ = [minPt.z,minPt.z+depth];

        // intervalX
        let intervalX;
        if (v.x < EPSILON && v.x>=0){
            if (boxX[0] <= e.x &&  e.x <= boxX[1]){
                intervalX = [-Infinity, Infinity];
            } else {
                return [];
            }
        } else {
            let aE = (boxX[0] - e.x)/v.x;
            let bE = (boxX[1] - e.x)/v.x;
            intervalX = [Math.min(aE,bE),Math.max(aE,bE)];
        }
        
        // intervalY
        let intervalY;
        if (v.y < EPSILON && v.y>=0){
            
            if (boxY[0] <= e.y &&  e.y <= boxY[1]){
                intervalY = [-Infinity, Infinity];
            } else {
                return [];
            }
        } else {
            let aE = (boxY[0] - e.y)/v.y;
            let bE = (boxY[1] - e.y)/v.y;
            intervalY = [Math.min(aE,bE),Math.max(aE,bE)];
        }
        
        // intervalZ
        let intervalZ;
        if (v.z < EPSILON && v.z>=0){
            if (boxZ[0] <= e.z && e.z <= boxZ[1]){
                intervalZ = [-Infinity, Infinity];
            } else {
                return [];
            }
        } else {
            let aE = (boxZ[0] - e.z)/v.z;
            let bE = (boxZ[1] - e.z)/v.z;
            intervalZ = [Math.min(aE,bE),Math.max(aE,bE)];
        }

        //check for overlap
        let allRight = this.overlap(intervalX, intervalY,intervalZ);
        if (allRight === undefined || (allRight[0]==allRight[1]) || Number.isNaN(allRight[0])){
            return [];
        }
        
        let t = allRight;
       
        const ret = [];
        for (const t of allRight){
            const pt = this.tToPt(t);
            let normal;
            if (Math.abs(pt.x - minPt.x) <= EPSILON) {
                normal = new Vector3(-1, 0, 0);
            } else if (Math.abs(pt.x - (minPt.x + g.v3_dim.x)) <= EPSILON) { 
                normal = new Vector3(1, 0, 0);
            } else if (Math.abs(pt.y - minPt.y) <= EPSILON) {
                normal = new Vector3(0, -1, 0);
            } else if (Math.abs(pt.y - (minPt.y + g.v3_dim.y)) <= EPSILON) { 
                normal = new Vector3(0, 1, 0);
            } else if (Math.abs(pt.z - minPt.z) <= EPSILON) {
                normal = new Vector3(0, 0, -1);
            } else if (Math.abs(pt.z - (minPt.z + g.v3_dim.z)) <= EPSILON) { 
                normal = new Vector3(0, 0, 1);
            } else {
                normal = vectorDifference(pt, minPt); 
            }
            normal.normalize();
            
            ret.push(new HitRecord(this,t,pt,g,normal));
        }
        return ret;
    }
    
    overlap(x,y,z){

        let left = Math.max(x[0],y[0],z[0]);
        let right = Math.min(x[1],y[1],z[1]);

        if (right > left){
            return [left,right];
        } else {
            return [0,0];
        }
    }
    
    
    hitCylinder(g){
        const center = g.v3_center;
        const height = g.f_height;
        const radius = g.f_radius;
        
        const testSphere = {
            v3_center: new Vector3(center.x,0,center.z),
            f_radius: radius,
        }
        const testRay = new Ray(new Vector3(this.start.x, 0, this.start.z), new Vector3(this.dir.x, 0, this.dir.z));
        const hits = testRay.hitSphere(testSphere);
        const ret = [];
        for (const h of hits){
            const pt = this.tToPt(h.t);
            if (center.y - height/2 < pt.y && pt.y < center.y + height/2){
                h.struckGeometry = g;
                h.pt = pt;
                h.ray = this;
                ret.push(h);
            }
        }
        
        for (const m of [1,-1]){
            const capCenter = new Vector3(center.x, center.y + m*height/2, center.z);
            const normal = new Vector3(0,m,0);
            const rhs = normal.dotProduct(capCenter);
            const t = (rhs - this.start.dotProduct(normal))/this.dir.dotProduct(normal);
            const hit = new HitRecord(this,t,this.tToPt(t),g,normal);
            const ptToAxis = vectorDifference(capCenter, hit.pt);
            if (ptToAxis.norm() < radius){
                ret.push(hit);
            }
        }
        return ret;
    }
    
    hitTriangle(_g){
        const pt0 = _g.v3_pt0; //a
        const pt1 = _g.v3_pt1; //b
        const pt2 = _g.v3_pt2; //c
        
        const start = this.start;
        const v = this.dir;
        
        let a = pt1.x - pt0.x; //bx-ax
        let b = pt1.y - pt0.y; //by-ay
        let c = pt1.z - pt0.z; //bz-az
        let d = pt2.x - pt0.x; //cx-ax
        let e = pt2.y - pt0.y; //cx-ay
        let f = pt2.z - pt0.z; //cx-az

        let g = -v.x;
        let h = -v.y; 
        let i = -v.z; 
        let j = start.x - pt0.x; //ex - ax
        let k = start.y - pt0.y; //ey - ay
        let l = start.z - pt0.z; //ez - az
          
        // solve M, t, lambda, beta
        const EIHF = e*i - h*f;
        const GFDI = g*f - d*i;
        const DHEG = d*h -e*g;
        const AKJB = a*k - j*b;
        const JCAL = j*c - a*l;
        const BLCK = b*l - k*c;
        
        const M = a*EIHF + b*GFDI + c*DHEG;        
        const t = -1*((f*AKJB + e*JCAL + d*BLCK)/M);

        const lambda = (i*AKJB + h*JCAL + g*BLCK)/M;
        if (lambda < 0 || lambda > 1){
            return [];
        }
        
        const beta = (j*EIHF + k*GFDI + l*DHEG)/M;
        if (beta < 0 || beta > 1 - lambda){
            return [];
        }
        
        // compute normal vector and return hit record
        const unit1 = vectorDifference(pt1, pt0);
        const unit2 = vectorDifference(pt2, pt0);
        const normal = unit1.crossProduct(unit2);
        normal.normalize();
        const ret = new HitRecord(this,t,this.tToPt(t),_g,normal);
        return [ret];        
    }
    
    hitsOfNode(AABBNode){
        // return list of hits between ray and objects contained at or below the given AABB-tree node

        let hit = this.hitBox(AABBNode.box);
        //if no bounding box
        if (hit.length == 0) {
            return [];
        }
        
        // if node is leaf
        if (AABBNode.geometries !== undefined) {
            return this.allHits(AABBNode.geometries);
        }
        const leftHits = this.hitsOfNode(AABBNode.left);
        const rightHits = this.hitsOfNode(AABBNode.right);
        return leftHits.concat(rightHits);    
    }    
}

class AABBNode {
    constructor(geometries) {
        // calls boxAround2(), boxAround(), and constructor in a loop
        
        // 1. bounding box
        let b = this.boxAround(geometries[0]);
        for (let i=1; i<geometries.length; i++){
            b = this.boxAround2(b,geometries[i]);    
        }
        this.box = b; //this box contains all the geometries within it
        
        // 2. if less than three object this is a leaf.
        if (geometries.length <= 3) {
            this.left = null;
            this.right = null;
            this.geometries = geometries;
            return;
        }
        
        // 3. If greater than three object we cut and split into left and right node.
        let maxDim = Math.max(this.box.v3_dim.x,this.box.v3_dim.y, this.box.v3_dim.z);
        if (maxDim === this.box.v3_dim.x){
            geometries.sort((a,b) => this.boxAround(a).v3_minPt.x - this.boxAround(b).v3_minPt.x);    
        } else if (maxDim === this.box.v3_dim.y){
            geometries.sort((a,b) => this.boxAround(a).v3_minPt.y - this.boxAround(b).v3_minPt.y);
        } else {
            geometries.sort((a,b) => this.boxAround(a).v3_minPt.z - this.boxAround(b).v3_minPt.z);
        }

        const len = geometries.length;
        const left = geometries.slice(0,len/2);
        const right = geometries.slice(len/2,len);
        
        this.left = new AABBNode(left);
        this.right = new AABBNode(right);
    }
    
    //return one geometry
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
    
    //two geometries
    boxAround2(g1,g2){       
        // create box for each object
        let box1 = this.boxAround(g1);
        let box2 = this.boxAround(g2);
        
        //calculate min and max
        let minX = Math.min(box1.v3_minPt.x, box2.v3_minPt.x);
        let minY = Math.min(box1.v3_minPt.y, box2.v3_minPt.y);
        let minZ = Math.min(box1.v3_minPt.z, box2.v3_minPt.z);
        let maxX = Math.max(box1.v3_minPt.x + box1.v3_dim.x, box2.v3_minPt.x + box2.v3_dim.x);
        let maxY = Math.max(box1.v3_minPt.y + box1.v3_dim.y, box2.v3_minPt.y + box2.v3_dim.y);
        let maxZ = Math.max(box1.v3_minPt.z + g1.v3_dim.z, box2.v3_minPt.z + box2.v3_dim.z);
        
        const ret = {
            s_type: 'box',
            v3_minPt: new Vector3(minX, minY, minZ),
            v3_dim: new Vector3(maxX-minX, maxY-minY, maxZ-minZ),
        };
        return ret;
        
    }
    
    boxPoints(pts) {
        /*
        Helper function: return box containing the given points
        */
        let [minX, minY, minZ] = [Infinity, Infinity, Infinity];
        let [maxX, maxY, maxZ] = [-Infinity, -Infinity, -Infinity];

        //check for bounds
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
        return this.boxPoints([
            g.v3_pt0,
            g.v3_pt1,
            g.v3_pt2, 
            vectorSum(g.v3_pt2, vectorDifference(g.v3_pt0, g.v3_pt1))
        ]);
    }
    
    boxCylinder(g) {
        const ret = {
            s_type: 'box',
            v3_minPt: vectorSum(g.v3_center, new Vector3(-g.f_radius, -g.f_height/2, -g.f_radius)),
            v3_dim: new Vector3(2*g.f_radius, g.f_height, 2*g.f_radius),
        };
        return ret;
    }

    
    boxSphere(g){
        const ret = {
            s_type: 'box',
            v3_minPt: vectorSum(g.v3_center, new Vector3(-g.f_radius, -g.f_radius, -g.f_radius)),
            v3_dim: new Vector3(2*g.f_radius, 2*g.f_radius, 2*g.f_radius),
        };
        return ret;    
    }
    
    boxTriangle(g){
        return this.boxPoints([
            g.v3_pt0,
            g.v3_pt1,
            g.v3_pt2,
            vectorSum(g.v3_pt2, vectorDifference(g.v3_pt0, g.v3_pt1))
        ]);    
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