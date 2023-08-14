import tensorflow as tf
import tensorflow.keras as keras
from tensorflow.keras import layers
from tensorflow.keras.callbacks import EarlyStopping 
from keras_bert import gelu

class CustomSchedule(tf.keras.optimizers.schedules.LearningRateSchedule):
  def __init__(self, d_model, warmup_steps=4000):
    super().__init__()

    self.d_model = d_model
    self.d_model = tf.cast(self.d_model, tf.float32)

    self.warmup_steps = warmup_steps

  def __call__(self, step):
    step = tf.cast(step, dtype=tf.float32)
    arg1 = tf.math.rsqrt(step)
    arg2 = step * (self.warmup_steps ** -1.5)

    return tf.math.rsqrt(self.d_model) * tf.math.minimum(arg1, arg2)

learning_rate = CustomSchedule(100)

opt = tf.keras.optimizers.Adam(learning_rate, beta_1=0.9, beta_2=0.98,
                                     epsilon=1e-9)

class TransformerBlock(layers.Layer):
    def __init__(self, embed_dim, num_heads, ff_dim, rate=0.1):
        super().__init__()
        self.att = layers.MultiHeadAttention(num_heads=num_heads, key_dim=embed_dim)
        self.ffn = keras.Sequential(
            [layers.Dense(ff_dim, activation="relu"), layers.Dense(embed_dim),]
        )
        self.layernorm1 = layers.LayerNormalization(epsilon=1e-6)
        self.layernorm2 = layers.LayerNormalization(epsilon=1e-6)
        self.dropout1 = layers.Dropout(rate)
        self.dropout2 = layers.Dropout(rate)

    def call(self, inputs, training):
        attn_output = self.att(inputs, inputs)
        attn_output = self.dropout1(attn_output, training=training)
        out1 = self.layernorm1(inputs + attn_output)
        ffn_output = self.ffn(out1)
        ffn_output = self.dropout2(ffn_output, training=training)
        return self.layernorm2(out1 + ffn_output)
    
class TokenAndPositionEmbedding(layers.Layer):
    def __init__(self, maxlen, vocab_size, embed_dim):
        super().__init__()
        self.token_emb = layers.Embedding(input_dim=vocab_size, output_dim=embed_dim)
        self.pos_emb = layers.Embedding(input_dim=maxlen, output_dim=embed_dim)

    def call(self, x):
        maxlen = tf.shape(x)[-1]
        print(maxlen)
        positions = tf.range(start=0, limit=maxlen, delta=1)
        positions = self.pos_emb(positions)
        x = self.token_emb(x)
        return x + positions


def def_model():
    embed_dim = 200  # Embedding size for each token
    num_heads = 2  # Number of attention heads
    ff_dim = 300
    maxlen = 300
    vocab_size = 4
    NUM_LAYERS = 2


    inputs = layers.Input(shape=(maxlen,))

    embedding_layer = TokenAndPositionEmbedding(maxlen, vocab_size, embed_dim)

    x = embedding_layer(inputs)
    conv_layer = layers.Conv1D(filters = 200,kernel_size=2, strides=1, padding="same", activation= gelu )
    #conv_layer2  = layers.Conv1D(filters = 100,kernel_size=5, strides=4, padding="same", activation= gelu)
    conv_layer3  = layers.Conv1D(filters = 100,kernel_size=7, strides=6, padding="same", activation= gelu)
    #conv_layer4  = layers.Conv1D(filters = 200,kernel_size=8, strides=7, padding="same", activation= gelu)
    #conv_layer5  = layers.Conv1D(filters = 200,kernel_size=8, strides=7, padding="same", activation= gelu)
    transformer_block = TransformerBlock(100, num_heads, ff_dim, rate=0.2)
    #transformer_block2 = TransformerBlock2(100, num_heads, ff_dim, rate=0.2)

    x = conv_layer(x)
    #x = conv_layer2(x)
    x = conv_layer3(x)
    #x = conv_layer4(x)
    #x = conv_layer5(x)
    for _ in range(NUM_LAYERS):
        x = TransformerBlock(100, num_heads, ff_dim, rate=0.2)(x)
    x = layers.GlobalAveragePooling1D()(x)
    x = layers.Dropout(0.6)(x)
    x = layers.Dense(20, activation="relu")(x)
    x = layers.Dropout(0.6)(x)

    outputs = layers.Dense(3, activation="softmax")(x)

    #categorical_crossentropy
    model = keras.Model(inputs=inputs, outputs=outputs)
    model.compile(optimizer=opt, loss='sparse_categorical_crossentropy', metrics=["sparse_categorical_accuracy"])
    model.summary()
    return model

