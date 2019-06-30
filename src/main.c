#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const unsigned short SceneGraphNode_AllocationStep = 10;

/* Types */
typedef float Vector4 __attribute((vector_size(16)));
typedef float Vector3 __attribute((vector_size(16)));
typedef float Vector2 __attribute((vector_size(8)));

typedef struct { Vector4 vector; } Quaternion;
typedef struct { Vector4 columns[4]; } Matrix;
typedef struct { Vector3 columns[3]; } Matrix3x3;

/* Operators */
inline const Vector3 Vector3Zero(void) {
    return (Vector3){ 0.f, 0.f, 0.f };
}

inline const Quaternion QuaternionIdentity(void) {
    return (Quaternion){ (Vector4){ 0, 0, 0, 1 } };
}

inline const Matrix MatrixIdentity(void) {
    return (Matrix) {
      (Vector4){ 1, 0, 0, 0 },
      (Vector4){ 0, 1, 0, 0 },
      (Vector4){ 0, 0, 1, 0 },
      (Vector4){ 0, 0, 0, 1 }
    };
}

#define SceneGraph_Radians(_degree) (_degree * M_PI / 180.f)
#define SceneGraph_Degree(_radians) (_radians * 180.0 / M_PI)

const float Vector3_Length(const Vector3 __x) {
  Vector3 _r = __x * __x;
  return _r[0] + _r[1] + _r[2];
}

const float Vector4_Length(const Vector4 __x) {
  Vector4 _r = __x * __x;
  return _r[0] + _r[1] + _r[2] + _r[3];
}

const Vector3 Vector3_Normalize(const Vector3 __x) { 
  return __x * (1 / sqrtf(Vector3_Length(__x)));
}

const Vector3 Vector3_Cross(const Vector3 __x, const Vector3 __y) { 
    Vector3 __xr = { __x[2], __x[0], __x[1] };
    Vector3 __yr = { __y[2], __y[0], __y[1] };
    Vector3 _r = __xr*__y - __x*__yr;
    return (Vector3){ _r[2], _r[0], _r[1] };
}

const float Vector3_Dot(const Vector3 __x, const Vector3 __y) {
    Vector3 _r = __x * __y;
  return _r[0] + _r[1] + _r[2];
}

const Matrix Matrix_Mul(const Matrix __x, const Matrix __y) { 
  Vector4 l11 = __x.columns[0][0] * __y.columns[0];
  Vector4 l12 = l11 + __x.columns[0][1] * __y.columns[1];
  Vector4 l13 = l12 + __x.columns[0][2] * __y.columns[2];
  Vector4 l14 = l13 + __x.columns[0][3] * __y.columns[3];

  Vector4 l21 = __x.columns[1][0] * __y.columns[0];
  Vector4 l22 = l21 + __x.columns[1][1] * __y.columns[1];
  Vector4 l23 = l22 + __x.columns[1][2] * __y.columns[2];
  Vector4 l24 = l23 + __x.columns[1][3] * __y.columns[3];

  Vector4 l31 = __x.columns[2][0] * __y.columns[0];
  Vector4 l32 = l31 + __x.columns[2][1] * __y.columns[1];
  Vector4 l33 = l32 + __x.columns[2][2] * __y.columns[2];
  Vector4 l34 = l33 + __x.columns[2][3] * __y.columns[3];

  Vector4 l41 = __x.columns[3][0] * __y.columns[0];
  Vector4 l42 = l41 + __x.columns[3][1] * __y.columns[1];
  Vector4 l43 = l42 + __x.columns[3][2] * __y.columns[2];
  Vector4 l44 = l43 + __x.columns[3][3] * __y.columns[3];

  return (Matrix) {
        (Vector4){ l14[0], l14[1], l14[2], l14[3] },
        (Vector4){ l24[0], l24[1], l24[2], l24[3] },
        (Vector4){ l34[0], l34[1], l34[2], l34[3] },
        (Vector4){ l44[0], l44[1], l44[2], l44[3] }
    }; 
}

const Vector3 Matrix3x3_Mul(const Vector3 __x, const Matrix3x3 __y) { 
  Vector3 l11 = __x[0] * __y.columns[0];
  Vector3 l12 = l11 + __x[1] * __y.columns[1];
  Vector3 l13 = l12 + __x[2] * __y.columns[2];
  return l13; 
}

Matrix Matrix_PerspectiveFOV(const float fovy, const float aspect, const float near, const float far) {
    const float angle  = SceneGraph_Radians(0.5f * fovy);
    const float yScale = 1.0f / tan(angle);
    const float xScale = yScale / aspect;
    const float zScale = far / (far - near);

    return (Matrix) {
        (Vector4){ xScale, 0.f, 0.f, 0.f },
        (Vector4){ 0.f, yScale, 0.f, 0.f },
        (Vector4){ 0.f, 0.f, zScale, 0.f },
        (Vector4){ 0.f, 0.f, -near * zScale, 0.f }
    };
}

Matrix Matrix_lookAt(const Vector3 eye, const Vector3 center, const Vector3 up)
{
    Vector3 zAxis = Vector3_Normalize(center - eye);
    Vector3 xAxis = Vector3_Normalize(Vector3_Cross(up, zAxis));
    Vector3 yAxis = Vector3_Cross(zAxis, xAxis);

    Vector4 P;
    Vector4 Q;
    Vector4 R;
    Vector4 S;

    P[0] = xAxis[0];
    P[1] = yAxis[0];
    P[2] = zAxis[0];
    P[3] = 0.0f;

    Q[0] = xAxis[1];
    Q[1] = yAxis[1];
    Q[2] = zAxis[1];
    Q[3] = 0.0f;

    R[0] = xAxis[2];
    R[1] = yAxis[2];
    R[2] = zAxis[2];
    R[3] = 0.0f;

    S[0] = -Vector3_Dot(xAxis, eye);
    S[1] = -Vector3_Dot(yAxis, eye);
    S[2] = -Vector3_Dot(zAxis, eye);
    S[3] =  1.0f;

    Matrix matrix = { P, Q, R, S };
    return matrix;
}

Matrix Matrix_Translate(const Vector3 vector) {
    Matrix matrix = MatrixIdentity();
    matrix.columns[3][0] = vector[0];
    matrix.columns[3][1] = vector[1];
    matrix.columns[3][2] = vector[2];
    return matrix;
}

const Vector3 Vector3_Invert(const Vector3 vector) {
    return (Vector3){ -vector[0], -vector[1], -vector[2] };
}

const Quaternion Quaternion_mul(Quaternion p, Quaternion q) {
    Vector4 _r = (p.vector[0] * __builtin_shufflevector(q.vector, -q.vector, 3, 6, 1, 4) +
                p.vector[1] * __builtin_shufflevector(q.vector, -q.vector, 2, 3, 4, 5)) +
                (p.vector[2] * __builtin_shufflevector(q.vector, -q.vector, 5, 0, 3, 6) +
                p.vector[3] * q.vector);
    return (Quaternion) { _r };
}

const Quaternion Quaternion_Inverse(Quaternion q) {
  return (Quaternion){ q.vector * (Vector4){ -1, -1, -1, 1 } * (1 / (Vector4_Length(q.vector))) };
}

const Matrix Quaternion_ToMatrix(Quaternion q) {
    Vector4 v = q.vector;
    Matrix r = {
    .columns[0] = { 1 - 2*(v[1]*v[1] + v[2]*v[2]),
                        2*(v[0]*v[1] + v[2]*v[3]),
                        2*(v[0]*v[2] - v[1]*v[3]), 0 },
    .columns[1] = {     2*(v[0]*v[1] - v[2]*v[3]),
                    1 - 2*(v[2]*v[2] + v[0]*v[0]),
                        2*(v[1]*v[2] + v[0]*v[3]), 0 },
    .columns[2] = {     2*(v[2]*v[0] + v[1]*v[3]),
                        2*(v[1]*v[2] - v[0]*v[3]),
                    1 - 2*(v[1]*v[1] + v[0]*v[0]), 0 },
    .columns[3] = { 0, 0, 0, 1 }
  };
  return r;
}

const Vector3 Vector3_MulMatrix(const Vector3 vector, const Matrix matrix) {
    Matrix3x3 m3x3 = { 
        (Vector3){ matrix.columns[0][0], matrix.columns[0][1], matrix.columns[0][2] },
        (Vector3){ matrix.columns[1][0], matrix.columns[1][1], matrix.columns[1][2] },
        (Vector3){ matrix.columns[2][0], matrix.columns[2][1], matrix.columns[2][2] },
    };
    return Matrix3x3_Mul(vector, m3x3);
}



/* SceneGraph */

struct SceneGraphNodes;
struct SceneGraphCamera;
struct SceneGraphNode;

// Структура, описывающая граф
typedef struct SceneGraph {
    /// Узлы графа в DoD формате
    struct SceneGraphNodes *nodes;

    /// Активная камера
    struct SceneGraphCamera *camera;
} SceneGraph;

typedef struct SceneGraphCameraUpdateCallback {
    struct SceneGraph *sceneGraph;
    void (*modifiedCallback)(struct SceneGraph *,
                             struct SceneGraphCamera *);
} SceneGraphCameraUpdateCallback;

typedef struct SceneGraphCamera {
    // View
    Vector3 position;
    Quaternion rotation;
    Vector3 target;
    Vector3 upVector;
    Matrix viewMatrix;

    // Projection
    float near;
    float far;
    float fov;
    float aspect;
    Matrix projectionMatrix;

    // Callbacks
    SceneGraphCameraUpdateCallback callback;
} SceneGraphCamera;

/**
 Уведомление о изменении матрицы модели
 */
typedef struct SceneGraphNodeUpdateCallback {
    struct SceneGraph *sceneGraph;
    void (*modifiedCallback)(struct SceneGraph *sceneGraph,
                             struct SceneGraphNode *node);
} SceneGraphNodeUpdateCallback;

/**
 Узлы графа в DoD формате
 */
typedef struct SceneGraphNodes {
    unsigned short capacity;
    unsigned short count;

    SceneGraphNodeUpdateCallback *callbacks;
    Vector3 *positions;
    Quaternion *rotations;
    Matrix *modelMatrixes;
    Matrix *MVPMatricies;
    unsigned char *isModified;
    struct SceneGraphNode **refs;
} SceneGraphNodes;

/**
 Указатель на узел графа в DoD формате
 */
typedef struct SceneGraphNode {
    SceneGraphNodes *nodes;
    unsigned short index;
} SceneGraphNode;


/* ---------- Camera ------------- */

static inline void _SceneGraphCamera_InvokeCallback(SceneGraphCamera *camera) {
    camera->callback.modifiedCallback(camera->callback.sceneGraph, camera);
}

static inline void _SceneGraphCamera_UpdateViewMatrix(SceneGraphCamera *camera) {
    camera->viewMatrix = Matrix_lookAt(camera->position,
                                         camera->position + camera->target,
                                         camera->upVector);
    _SceneGraphCamera_InvokeCallback(camera);
}

static inline void _SceneGraphCamera_UpdateProjectionMatrix(SceneGraphCamera *camera) {
    camera->projectionMatrix = Matrix_PerspectiveFOV(camera->fov,
                                                       camera->aspect,
                                                       camera->near,
                                                       camera->far);
    _SceneGraphCamera_InvokeCallback(camera);
}

#pragma mark - Camera

SceneGraphCamera *SceneGraphCamera_Create(SceneGraphCameraUpdateCallback callback) {
    SceneGraphCamera *camera = (SceneGraphCamera *)malloc(sizeof(SceneGraphCamera));
    camera->callback = callback;

    // View
    camera->position = Vector3Zero();
    camera->rotation = QuaternionIdentity();
    camera->target = (Vector3){ 0.f, 0.f, 1.f };
    camera->upVector = (Vector3){ 0.f, 1.f, 0.f };
    _SceneGraphCamera_UpdateViewMatrix(camera);

    // Projection
    camera->near = 0.1f;
    camera->far = 1000.f;
    camera->fov = 65.f;
    camera->aspect = 1.f;
    _SceneGraphCamera_UpdateProjectionMatrix(camera);

    return camera;
}

void SceneGraphCamera_Translate(SceneGraphCamera *camera, const Vector3 vector) {
    const Vector3 diff = vector - camera->position;
    const Matrix matrix = Matrix_Translate(Vector3_Invert(diff));
    camera->viewMatrix = Matrix_Mul(camera->viewMatrix, matrix);
    camera->position = camera->position + diff;
    _SceneGraphCamera_InvokeCallback(camera);
}

void SceneGraphCamera_Rotate(SceneGraphCamera *camera, const Quaternion quaternion) {
    const Quaternion diff = Quaternion_mul(quaternion, Quaternion_Inverse(camera->rotation));
    const Matrix matrix = Quaternion_ToMatrix(Quaternion_Inverse(diff));
    camera->viewMatrix = Matrix_Mul(matrix, camera->viewMatrix);
    camera->upVector = Vector3_MulMatrix(camera->upVector, matrix);
    camera->target = Vector3_MulMatrix(camera->target, matrix);
    camera->rotation = quaternion;
    _SceneGraphCamera_InvokeCallback(camera);
}

void SceneGraphCamera_SetAspect(SceneGraphCamera *camera, const float aspect) {
    camera->aspect = aspect;
    _SceneGraphCamera_UpdateProjectionMatrix(camera);
}

/* ------------- Node --------------- */

static inline void _SceneGraphNode_ApplyRotation(SceneGraphNode *node, const Quaternion quaternion) {
    const Matrix matrix = Quaternion_ToMatrix(quaternion);
    node->nodes->modelMatrixes[node->index] = Matrix_Mul(node->nodes->modelMatrixes[node->index], matrix);
}

static inline void _SceneGraphNode_ApplyTranslation(SceneGraphNode *node, const Vector3 vector) {
    const Matrix matrix = Matrix_Translate(vector);
    node->nodes->modelMatrixes[node->index] = Matrix_Mul(matrix, node->nodes->modelMatrixes[node->index]);
}

static inline void _SceneGraphNode_InvokeModifiedCallback(SceneGraphNode *node) {
    const SceneGraphNodeUpdateCallback *callBack = node->nodes->callbacks + node->index;
    callBack->modifiedCallback(callBack->sceneGraph, node);
}

#pragma mark - Nodes management

SceneGraphNodes *SceneGraphNode_CreateNodes() {
    SceneGraphNodes *nodes = (SceneGraphNodes *)malloc(sizeof(SceneGraphNodes));
    nodes->capacity = SceneGraphNode_AllocationStep;
    nodes->count = 0;
    nodes->positions = (Vector3 *)malloc(sizeof(Vector3) * nodes->capacity);
    nodes->rotations = (Quaternion *)malloc(sizeof(Quaternion) * nodes->capacity);
    nodes->callbacks = (SceneGraphNodeUpdateCallback *)malloc(sizeof(SceneGraphNodeUpdateCallback) * nodes->capacity);
    nodes->modelMatrixes = (Matrix *)malloc(sizeof(Matrix) * nodes->capacity);
    nodes->MVPMatricies = (Matrix *)malloc(sizeof(Matrix) * nodes->capacity);
    nodes->isModified = (unsigned char *)malloc(sizeof(unsigned char) * nodes->capacity);
    nodes->refs = (SceneGraphNode **)malloc(sizeof(SceneGraphNode *) * nodes->capacity);
    for (unsigned short i = 0; i < nodes->capacity; ++i) {
        nodes->refs[i] = (SceneGraphNode *)malloc(sizeof(SceneGraphNode));
        nodes->refs[i]->nodes = nodes;
        nodes->refs[i]->index = i;
    }
    return nodes;
}

SceneGraphNode *SceneGraphNode_Add(SceneGraphNodes *nodes, SceneGraphNodeUpdateCallback callback) {
    if (nodes->count >= nodes->capacity) {
        const unsigned short newCapacity = nodes->capacity + SceneGraphNode_AllocationStep;

        nodes->positions = (Vector3 *)realloc(nodes->positions, sizeof(Vector3) * newCapacity);
        nodes->rotations = (Quaternion *)realloc(nodes->rotations, sizeof(Quaternion) * newCapacity);
        nodes->callbacks = (SceneGraphNodeUpdateCallback *)realloc(nodes->callbacks, sizeof(SceneGraphNodeUpdateCallback) * newCapacity);
        nodes->modelMatrixes = (Matrix *)realloc(nodes->modelMatrixes, sizeof(Matrix) * newCapacity);
        nodes->MVPMatricies = (Matrix *)realloc(nodes->MVPMatricies, sizeof(Matrix) * newCapacity);
        nodes->isModified = (unsigned char *)realloc(nodes->isModified, sizeof(unsigned char) * newCapacity);
        nodes->refs = (SceneGraphNode **)realloc(nodes->refs, sizeof(SceneGraphNode *) * newCapacity);
        for (unsigned short i = nodes->capacity; i < newCapacity; ++i) {
            nodes->refs[i] = (SceneGraphNode *)malloc(sizeof(SceneGraphNode));
            nodes->refs[i]->nodes = nodes;
            nodes->refs[i]->index = i;
        }
        nodes->capacity = newCapacity;
    }

    SceneGraphNode *node = nodes->refs[nodes->count];
    nodes->callbacks[node->index] = callback;
    nodes->modelMatrixes[node->index] = MatrixIdentity();
    nodes->positions[node->index] = Vector3Zero();
    nodes->rotations[node->index] = QuaternionIdentity();
    nodes->isModified[node->index] = 0;
    nodes->count += 1;
    _SceneGraphNode_InvokeModifiedCallback(node);

    return node;
}

void SceneGraphNode_Delete(SceneGraphNodes *nodes, SceneGraphNode *node) {
    for (unsigned short i = node->index; i < nodes->count - 1; ++i) {
        nodes->positions[i] = nodes->positions[i + 1];
        nodes->rotations[i] = nodes->rotations[i + 1];
        nodes->callbacks[i] = nodes->callbacks[i + 1];
        nodes->modelMatrixes[i] = nodes->modelMatrixes[i + 1];
        nodes->MVPMatricies[i] = nodes->MVPMatricies[i + 1];
        nodes->isModified[i] = nodes->isModified[i + 1];
        SceneGraphNode *tmp = nodes->refs[i];
        nodes->refs[i] = nodes->refs[i + 1];
        nodes->refs[i + 1] = tmp;
        nodes->refs[i]->index = i;
    }
    nodes->count -= 1;
}

void SceneGraphNode_Translate(SceneGraphNode *node, const Vector3 vector) {
    _SceneGraphNode_ApplyTranslation(node, vector - node->nodes->positions[node->index]);
    node->nodes->positions[node->index] = vector;
    _SceneGraphNode_InvokeModifiedCallback(node);
}

void SceneGraphNode_Rotate(SceneGraphNode *node, const Quaternion quaternion) {
    const Quaternion diff = Quaternion_mul(quaternion, Quaternion_Inverse(node->nodes->rotations[node->index]));
    _SceneGraphNode_ApplyRotation(node, diff);
    node->nodes->rotations[node->index] = quaternion;
    _SceneGraphNode_InvokeModifiedCallback(node);
}

void SceneGraphNode_SetMVP(SceneGraphNode *node, const Matrix matrix) {
    node->nodes->MVPMatricies[node->index] = matrix;
}

const Matrix *SceneGraphNode_GetMVP(SceneGraphNode *node) {
    return node->nodes->MVPMatricies + node->index;
}

/* -------- SceneGraph --------- */
void _SceneGraph_UpdateNode(SceneGraphNode *node, SceneGraphCamera *camera) {
    node->nodes->MVPMatricies[node->index] = Matrix_Mul(Matrix_Mul(camera->projectionMatrix,
                                                                       camera->viewMatrix),
                                                          node->nodes->modelMatrixes[node->index]);
    node->nodes->isModified[node->index] = 0;
}

void _SceneGraph_CameraUpdateCallback(SceneGraph *sceneGraph, SceneGraphCamera *camera) {
    for (unsigned short i = 0; i < sceneGraph->nodes->count; ++i) {
        sceneGraph->nodes->isModified[i] = 1;
    }
}

void _SceneGraph_NodeUpdateCallback(SceneGraph *sceneGraph, SceneGraphNode *node) {
    node->nodes->isModified[node->index] = 1;
}

inline void _SceneGraph_Initialize(SceneGraph *sceneGraph) {
    sceneGraph->nodes = SceneGraphNode_CreateNodes();
    sceneGraph->camera = SceneGraphCamera_Create((SceneGraphCameraUpdateCallback) {
        .sceneGraph = sceneGraph,
        .modifiedCallback = _SceneGraph_CameraUpdateCallback
    });
}

SceneGraph *SceneGraph_Create() {
    SceneGraph *graph = (SceneGraph *)malloc(sizeof(SceneGraph));
    _SceneGraph_Initialize(graph);
    return graph;
}

void SceneGraph_Update(SceneGraph *sceneGraph) {
    for (unsigned short i = 0; i < sceneGraph->nodes->count; ++i) {
        if (sceneGraph->nodes->isModified[i] == 1) {
            _SceneGraph_UpdateNode(sceneGraph->nodes->refs[i], sceneGraph->camera);
        }
    }
}

SceneGraphNode *SceneGraph_CreateNode(SceneGraph *sceneGraph) {
    return SceneGraphNode_Add(sceneGraph->nodes, (SceneGraphNodeUpdateCallback) {
        .sceneGraph = sceneGraph,
        .modifiedCallback = _SceneGraph_NodeUpdateCallback
    });
}

void SceneGraph_DeleteNode(SceneGraph *sceneGraph, SceneGraphNode *node) {
    SceneGraphNode_Delete(sceneGraph->nodes, node);
}

SceneGraphCamera *SceneGraph_GetCamera(SceneGraph *sceneGraph) {
    return sceneGraph->camera;
}


/* Renderer */




/* WASM */
#define WASM_EXPORT __attribute__((visibility("default")))

/* External function that is implemented in JavaScript. */
extern void ClickAds();

extern void PrintMatrixLine(float x, float y, float z, float w);

/* Main */
WASM_EXPORT
SceneGraph *CreateScene() {
  SceneGraph *graph = SceneGraph_Create();
  SceneGraphCamera_Translate(graph->camera, (Vector3){ 5.f, 5.f, 5.f });

  int i, j;
  for (i = 0; i < 5; ++i) {
    for (j = 0; j < 5; ++j) {
      SceneGraphNode *node = SceneGraph_CreateNode(graph);
      SceneGraphNode_Translate(node, (Vector3){ i * 5.f, j * 5.f, 0.f });
    }
  }

  return graph;
}

WASM_EXPORT
void UpdateScene(SceneGraph *graph) {
  SceneGraph_Update(graph);
}

WASM_EXPORT
void DrawScene(SceneGraph *graph, void *context) {
    
}
