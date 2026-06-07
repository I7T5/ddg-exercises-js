"use strict";

/**
 * @module Projects
 */
class SimplicialComplexOperators {

        /** This class implements various operators (e.g. boundary, star, link) on a mesh.
         * @constructor module:Projects.SimplicialComplexOperators
         * @param {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
         * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
         */
        constructor(mesh) {
                this.mesh = mesh;
                this.assignElementIndices(this.mesh);

                this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
                this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
        }

        /** Assigns indices to the input mesh's vertices, edges, and faces
         * @method module:Projects.SimplicialComplexOperators#assignElementIndices
         * @param {module:Core.Mesh} mesh The input mesh which we index.
         */
        assignElementIndices(mesh) {
                // TODO
                indexElements(mesh.vertices); 
        }

        /** Returns the vertex-edge adjacency matrix of the given mesh.
         * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The vertex-edge adjacency matrix of the given mesh.
         */
        buildVertexEdgeAdjacencyMatrix(mesh) {
                // TODO
                let T = new Triplet(mesh.edges.length, mesh.vertices.length); 
                
                for (const e of mesh.edges) {
                        const r = e.index; 
                        const c1 = e.halfedge.vertex.index; 
                        const c2 = e.halfedge.twin.vertex.index; 
                        T.addEntry(1, r, c1); 
                        T.addEntry(1, r, c2); 
                }

                return SparseMatrix.fromTriplet(T); 
        }

        /** Returns the edge-face adjacency matrix.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeFaceAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The edge-face adjacency matrix of the given mesh.
         */
        buildEdgeFaceAdjacencyMatrix(mesh) {
                // TODO
                let T = new Triplet(mesh.faces.length, mesh.edges.length); 

                for (const f of mesh.faces) {
                        const r = f.index; 
                        for (let e of f.adjacentEdges()) {
                                const c = e.index; 
                                T.addEntry(1, r, c); 
                        }
                }

                return SparseMatrix.fromTriplet(T); 
        }

        /** Returns a column vector representing the vertices of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildVertexVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |V| entries. The ith entry is 1 if
         *  vertex i is in the given subset and 0 otherwise
         */
        buildVertexVector(subset) {
                // TODO
                let A = DenseMatrix.zeros(this.mesh.vertices.length); 
                for (const v of subset.vertices) {
                        A.set(1, v); 
                }
                return A; 
        }

        /** Returns a column vector representing the edges of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |E| entries. The ith entry is 1 if
         *  edge i is in the given subset and 0 otherwise
         */
        buildEdgeVector(subset) {
                // TODO
                let A = DenseMatrix.zeros(this.mesh.edges.length); 
                for (const e of subset.edges) {
                        A.set(1, e); 
                }
                return A; 
        }

        /** Returns a column vector representing the faces of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildFaceVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |F| entries. The ith entry is 1 if
         *  face i is in the given subset and 0 otherwise
         */
        buildFaceVector(subset) {
                // TODO
                let A = DenseMatrix.zeros(this.mesh.faces.length); 
                for (const f of subset.faces) {
                        A.set(1, f); 
                }
                return A; 
        }

        /** Returns the star of a subset.
         * 
         * The star St(S) is the collection of all simplices in K that contain any simplex in S.
         * 
         * @method module:Projects.SimplicialComplexOperators#star
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The star of the given subset.
         */
        star(subset) {
                // TODO
                let st = MeshSubset.deepCopy(subset);

                for (const v of this.mesh.vertices) {
                        if (st.vertices.has(v.index)) {
                                for (const e of v.adjacentEdges()) {
                                        st.addEdge(e.index); 
                                }
                        }
                }

                for (const e of this.mesh.edges) {
                        if (st.edges.has(e.index)) {
                                st.addFace(e.halfedge.face.index); 
                                st.addFace(e.halfedge.twin.face.index); 
                        }
                }

                // The star of a face is the emptyset. 

                return st; 
        }

        /** Returns the closure of a subset.
         * 
         * The closure Cl(S) is the smallest (i.e., fewest elements) subcomplex of K that contains S. 
         * 
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
                // TODO
                let cl = MeshSubset.deepCopy(subset);

                for (const f of this.mesh.faces) {
                        if (cl.faces.has(f.index)) {
                                for (const e of f.adjacentEdges()) {
                                        cl.addEdge(e.index); 
                                }  
                        }
                }

                for (const e of this.mesh.edges) {
                        if (cl.edges.has(e.index)) {
                                cl.addVertex(e.halfedge.vertex.index); 
                                cl.addVertex(e.halfedge.twin.vertex.index); 
                        }
                }

                return cl; 
        }

        /** Returns the link of a subset.
         * 
         * The link Lk(S) is equal to Cl(St(S)) \ St(Cl(S)).
         * 
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {
                // TODO
                const vs = this.closure(this.star(subset)).vertices.difference(this.star(this.closure(subset)).vertices); 
                const es = this.closure(this.star(subset)).edges.difference(this.star(this.closure(subset)).edges); 
                const fs = this.closure(this.star(subset)).faces.difference(this.star(this.closure(subset)).faces); 
                
                return new MeshSubset(vs, es, fs); 
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {
                // TODO
                // A subset is a simplicial complex if it's equal to its closure. 
                // Every vertex contains itself. 
                // Every edge contains its adjacent vertices. 
                // Every face contains its adjacent edges. 
                return subset.equals(this.closure(subset)); 
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {
                // TODO
                const emptyset = new Set(); 
                if (subset.faces.size > 0) {
                        let faceSubset = new MeshSubset(emptyset, emptyset, subset.faces); 

                        // Pure if subset === Cl(justTheFaces)
                        if (this.closure(faceSubset).equals(subset)) return 2; 
                        return -1; 
                } 

                if (subset.edges.size > 0) {
                        let edgeSubset = new MeshSubset(emptyset, subset.edges, emptyset); 

                        if (this.closure(edgeSubset).equals(subset)) return 1; 
                        return -1; 
                }

                return 0; 
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
                // TODO
                let bd = new MeshSubset(); 
                switch (this.isPureComplex(subset)) {
                        case 0: 
                                break; 
                        case 1: 
                                // throw away duplicate vertices of edge.adjVertices
                                for (const e of this.mesh.edges) {
                                        if (subset.edges.has(e.index)) {
                                                const vs = [e.halfedge.vertex, e.halfedge.twin.vertex]; 
                                                for (const v of vs) {
                                                        if (!bd.vertices.has(v.index)) {
                                                                bd.addVertex(v.index); 
                                                        } else {
                                                                bd.deleteVertex(v.index); 
                                                        }
                                                }
                                        }
                                }
                                break; 
                        case 2: 
                                // throw away duplicate edges of face.adjEdges
                                for (const f of this.mesh.faces) {
                                        if (subset.faces.has(f.index)) {
                                                for (const e of f.adjacentEdges()) {
                                                        if (!bd.edges.has(e.index)) {
                                                                bd.addEdge(e.index); 
                                                        } else {
                                                                bd.deleteEdge(e.index); 
                                                        }
                                                }
                                        }
                                }

                                // vertices are the adjacent vertices of remaining edges
                                for (const e of this.mesh.edges) {
                                        if (subset.edges.has(e.index)) {
                                                bd.addVertex(e.halfedge.vertex.index); 
                                                bd.addVertex(e.halfedge.twin.vertex.index); 
                                        }
                                }
                                break; 
                        
                }

                return bd; 
        }
}
